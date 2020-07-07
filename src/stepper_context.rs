use acc_vector::AccSegment;
use acc_vector::AccVector;
use integer_curve::PathLimits;
use integer_curve;
//use bezier;
use num::Integer;
use std::u64;
//use std::convert::TryFrom;
use curve_approx::CurveInfo;
use coords::{Point, Vector};
use curves;
use std::error::Error;


#[derive(Debug)]
pub enum CurveSegment
{
    GoTo(Point), // x,y
    GoToRel(Point), // x,y
    LineTo(Point), // x,y
    LineToRel(Vector), // x,y
    // Control points are relative to the closest end point
    CurveTo(Point, Vector, Vector), // p2, c1, c2
    CurveToRel(Vector, Vector, Vector), // p2, c1, c2
    Arc(f64, f64, f64, f64, f64) // rx,ry, start, end, rotation
}


#[derive(Debug, Clone)]
pub enum Command
{
    Acc(u8, i16), // (channel, acceleration)
    Weight(i32) // (Weight value)
}

#[derive(Debug, Clone)]
pub struct StepperEvent
{
    pub ticks :u32, // Ticks since previous event
    pub cmd: Command
}

const N_CHANNELS : usize = 2;
const X_INDEX : usize = 0;
const Y_INDEX : usize = 1;

/// State for a single channel (stepper)
#[derive(Copy, Clone, Debug)]
pub struct ChannelContext
{
    pos : i64,  // Position in steps
    v : i32, // steps/(2*tick)
    ticks : u64 // Ticks so far in this channel
}

/// State for curve generation
#[derive(Debug)]
pub struct StepperContext
{
    /// State for each channel
    channels: [ChannelContext; N_CHANNELS],
    /// Maximum acceleration for each channel
    step_a_max: [i32; N_CHANNELS],
    /// Maximum speed for each channel
    step_v_max: [i32; N_CHANNELS],

    events: Vec<StepperEvent>,
    event_ticks: u64, // Time of last event

    move_weight: i32,
    curve_weight: i32,

    /// Conversion factor for coordinates in curves to steps
    step_scale: [f64; N_CHANNELS],
    
    /// Conversion factor for velocity in curves (units/tick) to steps (steps/(2*tick))
    v_scale: [f64; N_CHANNELS],

    
    // Minimum cosine value for an angle between endpoints of two curves to
    // allow connecting then
    min_cos_connect: f64
}


fn curve_segment_to_info(seg: &CurveSegment, current_pos: &mut Point)
                         -> Option<Box<dyn CurveInfo>>
{
    Some(
        match seg {
            CurveSegment::LineTo(p2) =>  {
                let rel = *p2 - *current_pos;
                *current_pos = *p2;
                Box::new(curves::line::Line::new(rel))
            },
            CurveSegment::LineToRel(p2) => {
                *current_pos += *p2;
                Box::new(curves::line::Line::new(*p2))
            },
            CurveSegment::CurveTo(p2, c1, c2) =>  {
                let rel = *p2 - *current_pos;
                *current_pos = *p2;
                Box::new(curves::bezier::Bezier::new(*c1, rel + *c2, rel))
            },
            CurveSegment::CurveToRel(p2, c1 ,c2) => {
                *current_pos += *p2;
                Box::new(curves::bezier::Bezier::new(*c1, *p2 + *c2, *p2))
            },
            CurveSegment::Arc(rx,ry, start, end, _rot)  => {
                // Only circles are supported
                if *rx != *ry {return None;}
                let circle = curves::circle_segment::CircleSegment::new(*rx, *start, *end);
                let (end, _) = circle.value(circle.length());
                *current_pos += end;
                Box::new(circle)
            }
            _ => return None
        })
}

fn limit<T>(a: T, min: T, max: T) -> T
    where T: PartialOrd
{
    if a > max {
        max
    } else if a < min {
        min
    } else {
        a
    }
}

fn acc_max<T>(acc: &mut T, b: T)
    where T: PartialOrd
{
    if *acc < b {
        *acc = b;
    }
}

fn acc_max_err(acc: &mut [i64;N_CHANNELS], err: &[i64;N_CHANNELS])
{
    for d in 0..N_CHANNELS {
        acc_max(&mut acc[d], err[d]);
    }
}

// Allow up to 5 degrees angle between end points
const INITIAL_MIN_COS_CONNECT: f64 = 0.99619;

impl StepperContext {
    pub fn new(a_max: &[i32;N_CHANNELS], v_max: &[i32;N_CHANNELS],
               step_scale: &[f64;N_CHANNELS],  v_scale: &[f64;N_CHANNELS])
               -> StepperContext
    {
        StepperContext {channels: [ChannelContext {pos:0,v:0,ticks:0} ; 
                                   N_CHANNELS] ,
                        step_a_max: *a_max, step_v_max: *v_max,
                        events: Vec::<StepperEvent>::new(),
                        event_ticks : 0,
                        move_weight: 0,
                        curve_weight: 255,
                        step_scale: *step_scale,
                        v_scale: *v_scale,
                        min_cos_connect: INITIAL_MIN_COS_CONNECT
        }
    }

    pub fn add_acc_interval(&mut self, acc_x: i16, acc_y: i16, interval: u32)
    {
        let acc = [acc_x, acc_y];
        assert!(interval >= 1);
        for i in 0..N_CHANNELS {
            /*
            if self.channels[i].ticks == self.event_ticks {
                self.events.pop();
            }
             */
            self.events.push(StepperEvent {
                ticks: (self.channels[i].ticks-self.event_ticks) as u32,
                cmd: Command::Acc(i as u8, acc[i])
            });
            self.channels[i].pos += 
                2 * interval as i64 * self.channels[i].v as i64
                + interval as i64 * interval as i64 * acc[i] as i64;
            self.channels[i].v += interval as i32 * acc[i] as i32;
            self.event_ticks = self.channels[i].ticks;
            self.channels[i].ticks += interval as u64;
        }
        assert_eq!(self.channels[X_INDEX].ticks, self.channels[Y_INDEX].ticks);
    }

    pub fn add_weight_change(&mut self, weight: i32, when: u64)
    {
        assert!(when >= self.channels[X_INDEX].ticks);
        assert!(when >= self.channels[Y_INDEX].ticks);
        self.events.push(StepperEvent {
            ticks: (when-self.event_ticks) as u32,
            cmd: Command::Weight(weight)
        });
        self.event_ticks = when;
    }
    
    fn acc_segs_to_events(&mut self, acc: &[Vec<AccSegment>]) 
    {
        let mut p = Vec::new();
        for c in acc {
            p.push(c.iter().peekable());
        }
        loop {
            let mut min = u64::MAX;
            let mut min_index : Option<usize> = None;
            for i in 0..N_CHANNELS {
                if self.channels[i].ticks < min {
                    if p[i].peek().is_some() {
                        min_index = Some(i);
                        min = self.channels[i].ticks;
                    }
                }
            }
            if let Some(i) = min_index {
                if let Some(seg) = p[i].next() {
                    assert!(self.channels[i].ticks >= self.event_ticks);
                    self.events.push(StepperEvent {
                        ticks: (self.channels[i].ticks-self.event_ticks) as u32,
                        cmd: Command::Acc(i as u8, seg.acc)
                    });
                    self.event_ticks = self.channels[i].ticks;
                    self.channels[i].ticks += seg.interval as u64;
                }
            } else {
                break;
            }
        }
    }

    /*
    fn event_distance(&self, end: u64) ->([i64;N_CHANNELS], [i32; N_CHANNELS])
    {
        let mut a = [0i16;N_CHANNELS];
        let mut v = [0i32;N_CHANNELS];
        let mut s = [0i64;N_CHANNELS];
        let mut t_event = 0u64;
        for e in &self.events {
            t_event += e.ticks as u64;
            if t_event >= end {
                t_event -= e.ticks as u64;
                break;
            }
            for ch in 0..N_CHANNELS {
                let v_next = v[ch] + a[ch] as i32 * e.ticks as i32;
                s[ch] += (v_next + v[ch]) as i64 * e.ticks as i64;
                v[ch] = v_next;
            }
            match e.cmd {
                Command::Acc(ch,acc) => {
                    a[ch as usize] = acc;
                },
                _ => {}
            }
        }
        for ch in 0..N_CHANNELS {
            let dt = end - t_event;
            let v_next = v[ch] + a[ch] as i32 * dt as i32;
            s[ch] += (v_next + v[ch]) as i64 * dt as i64;
            v[ch] = v_next;
        }
        
        (s, v)
    }
   
     */
    fn to_steps(&self, coord: f64, channel: usize) -> i64
    {
        let s = (coord*self.step_scale[channel]).round();
        assert!(s < (i64::max_value() as f64)
                && s > (i64::min_value() as f64),
                "Step coordinate out of range");
        s as i64
    }
     
    fn to_step_v(&self, v: f64, channel: usize) -> i32
    {
        let s = (v * self.v_scale[channel]).round();
        assert!(s < (i32::max_value() as f64)
                && s > (i32::min_value() as f64),
                "Step coordinate out of range");
        s as i32
    }


    fn from_steps(&self, steps: i64, channel: usize) -> f64
    {
        (steps as f64) / self.step_scale[channel]
    }
    
    
    /// Goto a position given as steps, with given final velocity
    pub fn step_goto_speed(&mut self, x: i64, y: i64, vx: i32, vy:i32)
    {
        assert!(self.step_v_max[X_INDEX] >= vx);
        assert!(self.step_v_max[Y_INDEX] >= vy);
        self.events.push(StepperEvent {
            ticks: (self.channels[X_INDEX].ticks-self.event_ticks) as u32,
            cmd: Command::Weight(self.move_weight)
        });
        self.event_ticks = self.channels[X_INDEX].ticks;
        
        let limits = &[PathLimits {ds: x - self.channels[X_INDEX].pos,
                                   v0: self.channels[X_INDEX].v, vn: vx, 
                                   a: self.step_a_max[X_INDEX], 
                                   vmax: self.step_v_max[X_INDEX]},
                       PathLimits {ds: y - self.channels[Y_INDEX].pos,
                                   v0: self.channels[Y_INDEX].v, vn: vy, 
                                   a: self.step_a_max[Y_INDEX],
                                   vmax: self.step_v_max[Y_INDEX]}
        ];
        //println!("limits: {:?}", *limits);
        let seq = integer_curve::shortest_curve_sequences(limits);
        let v_end = [vx,vy];
        for i in 0..N_CHANNELS {
            //println!("{}: {:?}", i, seq[i]);
            let (s, v) = seq[i].acc_distance(self.channels[i].v);
            self.channels[i].pos += s;
            assert_eq!(v, v_end[i]);
        }
        self.channels[X_INDEX].v = vx;
        self.channels[Y_INDEX].v = vy;
        self.acc_segs_to_events(&seq);
        assert_eq!(self.channels[X_INDEX].ticks, self.channels[Y_INDEX].ticks);

        /*
        let (s,v) = self.event_distance(self.channels[X_INDEX].ticks);
        println!("s={:?}, v={:?}", s,v);
         */
    }
    
    pub fn step_goto(&mut self, x: i64, y: i64)
    {
        self.step_goto_speed(x,y,0,0)
    }

    
    pub fn goto_speed(&mut self, x: f64, y: f64, vx: f64, vy:f64)
    {
        self.step_goto_speed(self.to_steps(x, X_INDEX),
                             self.to_steps(y, Y_INDEX), 
                             self.to_step_v(vx, X_INDEX),
                             self.to_step_v(vy, Y_INDEX));
    }
    
     pub fn goto(&mut self, x: f64, y: f64)
    {
        self.step_goto_speed(self.to_steps(x, X_INDEX),
                             self.to_steps(y, Y_INDEX),
                             0,0)
    }

    pub fn speed(&mut self, v_x: f64, v_y: f64) {
        let v = [self.to_step_v(v_x,X_INDEX), self.to_step_v(v_y,Y_INDEX)];
        let mut acc = [Vec::<AccSegment>::new(), Vec::<AccSegment>::new()];
        for i in 0..N_CHANNELS {
            let (t,dt) = 
                (v[i] - self.channels[i].v).div_rem(&self.step_a_max[i]);
            let a = if v[i] >= self.channels[i].v {
                self.step_a_max[i]
            } else {
                -self.step_a_max[i]
            };
            acc[i].push(AccSegment {interval: t.abs() as u16, 
                                    acc: a as i16});
            self.channels[i].pos += 
                ((2*self.channels[i].v + self.step_a_max[i] * t) * t.abs())
                as i64;
            if dt != 0 {
                acc[i].push(AccSegment {interval: 1, acc: dt as i16});
                self.channels[i].pos += (self.channels[i].v 
                                         + self.step_a_max[i] * t + v[i])
                    as i64;
            }
            self.channels[i].v = v[i]; 
        }
        self.acc_segs_to_events(&acc);
    }

    
    pub fn approx_curve(&mut self, start: Point, curve: &dyn CurveInfo, v: f64)
        -> [i64; 2]
    {
        let len = curve.length();
        let steps = (len/v).round();
        let step_len = len / steps;
        let weigth = 0.5;
        let mut max_err = [0i64,0i64];
        
        let mut pos_err_acc = [0.0, 0.0];
        
        let a_max = [self.step_a_max[0] as f64, self.step_a_max[1] as f64];
        let v_max = [self.step_v_max[0], self.step_v_max[1]];
        for t in 1..=(steps as i32) {
            let p = (t as f64) * step_len;
            let (pos,dir) = curve.value(p);
            let mut acc = [0i16;N_CHANNELS];
            /*
            println!("p: {} pos: {} ({}, {}) ",p, pos+start,
                     (pos+start).x * self.step_scale[X_INDEX],
                     (pos+start).y * self.step_scale[Y_INDEX]);
             */
            let mut step_pos_int = [0i64, 0i64];
            for d in 0..N_CHANNELS {
                let step_v = dir[d] * v * self.v_scale[d];
                let step_pos = (pos[d]+start[d]) * self.step_scale[d];
                step_pos_int[d] = step_pos as i64;
                // Acceleration for correct position
                let a_pos = limit(
                    (step_pos - (self.channels[d].pos
                                + 2* self.channels[d].v as i64) as f64) / 2.0,
                    -a_max[d], a_max[d]);
                pos_err_acc[d] += a_pos * 0.1;

                // Clear accumulator if the error changed sign
                if (pos_err_acc[d] > 0.0 && a_pos < 0.0) 
                    || (pos_err_acc[d] < 0.0 && a_pos > 0.0)
                {
                    pos_err_acc[d] = 0.0;
                }
                
                // Acceleration for correct speed
                let a_v = limit(
                    step_v - self.channels[d].v as f64,
                    -a_max[d], a_max[d]);
                
                // println!("a_pos: {} a_v: {}", a_pos, a_v);
                acc[d] = (a_v * (1.0-weigth) 
                          + (a_pos + pos_err_acc[d]) * weigth).round() as i16;
                
                // Limit speed
                if acc[d] < 0 
                    && self.channels[d].v + (acc[d] as i32) < -v_max[d] 
                {
                    acc[d] = (-v_max[d] - self.channels[d].v) as i16;
                } else if acc[d] > 0 
                    && self.channels[d].v + (acc[d] as i32) > v_max[d] 
                {
                    acc[d] = (v_max[d] - self.channels[d].v) as i16;
                }
            }
            self.add_acc_interval(acc[X_INDEX], acc[Y_INDEX], 1);
            acc_max(&mut max_err[X_INDEX], 
                    (self.channels[X_INDEX].pos - step_pos_int[X_INDEX]).abs());
            acc_max(&mut max_err[Y_INDEX], 
                    (self.channels[Y_INDEX].pos - step_pos_int[Y_INDEX]).abs());
            //println!("channel: {:?} ", self.channels);
            
        }
        max_err
    }
    
   
    
    pub fn set_weight(&mut self, weight: i32)
    {
        self.curve_weight = weight;
    }

 

    fn draw_curve(&mut self, start: Point, curve: &mut curves::concat_curve::ConcatCurve, v: f64) -> Result<[i64;2], Box<dyn Error>>
    {
        if !curve.is_empty() {
            println!("Approximating curve: {} {:?}", start, curve);
            let (_, start_dir) = curve.value(0.0);
            let v_start = start_dir * v;
            self.add_weight_change(self.move_weight,
                                   self.channels[X_INDEX].ticks);
            self.goto_speed(start.x, start.y, v_start.x, v_start.y);
            self.add_weight_change(self.curve_weight,
                                   self.channels[X_INDEX].ticks);
            let max_err = self.approx_curve(start, curve, v);
            curve.clear();
            Ok(max_err)
        } else {
            Ok([0,0])
        }
    }

    
    pub fn draw_curves(&mut self, segs: &[CurveSegment], v: f64)
        -> Result<[i64;2], Box<dyn Error>>
    {
        let mut max_err = [0i64,0i64];
        let mut iter = segs.iter();
        let mut seg =
            match iter.next() {
                Some(s) => s,
                None => return Ok(max_err)
            };
        let mut current_pos = Point {
            x: self.from_steps(self.channels[X_INDEX].pos, X_INDEX),
            y: self.from_steps(self.channels[Y_INDEX].pos, Y_INDEX)
        };
        let mut start = current_pos;
        let mut prev_dir: Option<Vector> = None;
        let mut curves = curves::concat_curve::ConcatCurve::new();
        loop {
            let mut next_pos = current_pos;
            if let Some(info) = curve_segment_to_info(&seg, &mut next_pos) {
                let (_,start_dir) = info.value(0.0);
                if let Some(prev_dir) = prev_dir {
                    println!("Direction: {} -> {}", prev_dir, start_dir);
                    if start_dir.x * prev_dir.x + start_dir.y * prev_dir.y 
                        < self.min_cos_connect {
                            println!("Splitting curve");
                            let err = self.draw_curve(start, &mut curves, v)?;
                            acc_max_err(&mut max_err, &err);
                            start = current_pos;
                        }
                }
                let (_, end_dir) = info.value(info.length());
                current_pos = next_pos;
                prev_dir = Some(end_dir);
                curves.add(info);
            } else {
                match *seg {
                    CurveSegment::GoTo(p) => {
                        // Skip gotos to the current position
                        if self.to_steps(p.x, X_INDEX) 
                            != self.channels[X_INDEX].pos 
                            || self.to_steps(p.y, Y_INDEX)
                            != self.channels[Y_INDEX].pos {
                                let err=self.draw_curve(start, &mut curves, v)?;
                                acc_max_err(&mut max_err, &err);
                                current_pos = p;
                                start = current_pos;
                                prev_dir = None;
                            }
                    },
                    CurveSegment::GoToRel(rel) =>  {
                        if rel.x != 0.0 || rel.y != 0.0 {
                            let err = self.draw_curve(start, &mut curves, v)?;
                            acc_max_err(&mut max_err, &err);
                            current_pos += rel;
                            start = current_pos;
                            prev_dir = None;
                        }
                    },
                    _ => panic!("Unhandled CurveSegment type")
                };
            }
            
            let next_seg = iter.next();
            if let Some(s) = next_seg {
                seg = s;
            } else {
                break;
            }
            println!("Current pos: {}", current_pos);
        }
        let err = self.draw_curve(start, &mut curves, v)?;
        acc_max_err(&mut max_err, &err);
        self.add_weight_change(self.move_weight, self.channels[X_INDEX].ticks);
        Ok(max_err)
    }

    pub fn position(&self) -> (i64, i64)
    {
        (self.channels[X_INDEX].pos, self.channels[Y_INDEX].pos)
    }
    
    pub fn velocity(&self) -> (i32, i32)
    {
        (self.channels[X_INDEX].v, self.channels[Y_INDEX].v)
    }
    
    pub fn ticks(&self) -> u64
    {
        assert_eq!(self.channels[X_INDEX].ticks, self.channels[Y_INDEX].ticks);
        self.channels[X_INDEX].ticks
    }

    
    pub fn events<'a>(&'a mut self) -> &'a Vec<StepperEvent>
    {
        assert_eq!(self.channels[X_INDEX].ticks, self.channels[Y_INDEX].ticks);

        
        self.events.push(StepperEvent {
            ticks: (self.channels[X_INDEX].ticks-self.event_ticks) as u32,
            cmd: Command::Weight(self.move_weight)
        });
        self.event_ticks = self.channels[X_INDEX].ticks;
        
        for i in 0..N_CHANNELS {
            self.events.push(StepperEvent {
                ticks: (self.channels[i].ticks-self.event_ticks) as u32,
                cmd: Command::Acc(i as u8, 0)
            });
            self.event_ticks = self.channels[i].ticks;
        }
        &self.events
    }
}

#[test]
fn test_approx_curve_line()
{
    let mut ctxt = StepperContext::new(&[10, 10],
                                       &[10, 10],
                                       &[3.0, 7.0],
                                       &[1.5, 3.5]);

    let line = curves::line::Line::new(Point{x: 23.0, y: 45.0});
    ctxt.approx_curve(Point{x: 0.0, y:0.0}, &line, 1.0);

}

#[cfg(test)]
use std::f64::consts::PI;

#[test]
fn test_approx_curve_ellipse()
{
    let mut ctxt = StepperContext::new(&[10, 10],
                                       &[10, 10],
                                       &[3.0, 7.0],
                                       &[1.5, 3.5]);

    let circle = curves::circle_segment::CircleSegment::new(670.0, 0.0, PI/2.0);
    ctxt.approx_curve(Point{x: 0.0, y:0.0}, &circle, 30.0);

}
