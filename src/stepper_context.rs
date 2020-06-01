use acc_vector::AccSegment;
use acc_vector::AccVector;
use integer_curve::PathLimits;
use integer_curve;
//use bezier;
use num::Integer;
use std::u64;
//use std::convert::TryFrom;
use curve_approx::CurveInfo;
use std::f64::consts::PI;

#[derive(Debug)]
pub enum CurveSegment
{
    GoTo(f64, f64), // x,y
    GoToRel(f64, f64), // x,y
    LineTo(f64, f64), // x,y
    LineToRel(f64, f64), // x,y
    // Control points are relative to the closest end point
    CurveTo(f64, f64, f64, f64, f64, f64), // x,y, c1x, c1y, c2x, c2y
    CurveToRel(f64, f64, f64, f64, f64, f64), // x,y, c1x, c1y, c2x, c2y
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
    /// Conversion factor for velocity in curves (units/s) to steps (steps/(2*tick))
    v_scale: [f64; N_CHANNELS],
}

// Returns (vx, vy)
fn vel_dir(x:f64, y:f64, v:i32) -> (i32, i32)
{
    let dist = (x*x + y*y).sqrt();
    let scale = (v as f64) / dist;
    ((x * scale).round() as i32,
     (y * scale).round() as i32)
}

const ANGLE_SIN_LIMIT: f64 = 0.1; 
const SPEED_REL_LIMIT: f64 = 0.10; // Speeds within 10%


// Check if velocities matches in direction and magnitude
fn velocity_matches(v1x:f64, v1y:f64, v2x:f64, v2y:f64,
                    angle_sin_limit:f64, speed_rel_limit:f64) -> bool
{
    let v1 = (v1x*v1x + v1y*v1y).sqrt();
    let v2 = (v2x*v2x + v2y*v2y).sqrt();
    if v1 == 0.0 || v2 == 0.0 {
        return v1 == v2;
    }
    if if v1 > v2 { v2 / v1 } else { v1 / v2} < (1.0-speed_rel_limit) {
        println!("Speed mismatch");
        // Difference in magnitude too big
        false
    } else {
        if (v1x*v2y - v1y*v2x).abs() / (v1 * v2) < angle_sin_limit {
            println!("Angle mismatch: {}", (v1x*v2y - v1y*v2x));
        }
        (v1x*v2y - v1y*v2x).abs() / (v1 * v2) < angle_sin_limit
    }
}

struct Line
{
    x: f64,
    y: f64,
    dir: [f64;2]
}

impl Line
{
    pub fn new( x: f64, y: f64) -> Line
    {
        let l = (x*x + y*y).sqrt();
        Line{x,y, dir: [x/l, y/l]}
    }
}
impl CurveInfo for Line
{
    fn length(&self) -> f64
    {
        (self.x*self.x + self.y*self.y).sqrt()
    }
    fn value(&self, pos: f64) -> ([f64;2], [f64;2])
    {
        ([self.x*pos, self.y*pos], self.dir)
    }
}

struct Ellipse
{
    rx: f64,
    ry: f64,
    start: f64,
    arc: f64,
    rotation: f64,
    x0: f64,
    y0: f64
}
impl Ellipse
{
    pub fn new(rx: f64,
               ry: f64,
               start: f64,
               end: f64,
               rotation: f64) -> Ellipse
    {
        if rx != ry {
            panic!("Non circular ellipses not supported");
            // TODO
        }
        let mut arc = end - start;
        if arc <= 0.0 {
            arc += PI;
        }
        let (s, c) = start.sin_cos();
        let x0 = -c * rx;
        let y0 = -s * ry;
        Ellipse{rx, ry, start, arc, rotation, x0, y0}
    }
}

impl CurveInfo for Ellipse
{
    fn length(&self) -> f64
    {
        self.rx * self.arc
    }
    
    fn value(&self, pos: f64) -> ([f64;2], [f64;2])
    {
        let (s, c) = (self.start + self.arc * pos).sin_cos();
        ([self.x0 + c * self.rx, self.y0 + s * self.ry], [-s, c])
    }
}



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
                        v_scale: *v_scale
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
        assert!(when >= self.event_ticks);
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

    
    pub fn approx_curve(&mut self, curve: &dyn CurveInfo, v: f64)
    {
        let len = curve.length();
        let steps = (len/v).round();
        let step_len = 1.0 / steps;
        let s0 = [self.channels[X_INDEX].pos, self.channels[Y_INDEX].pos];
        for t in 1..(steps as i32 + 1) {
            let p = (t as f64) * step_len;
            let (pos,dir) = curve.value(p);
            let mut acc = [0i16;N_CHANNELS];
            println!("p: {} pos: {:?} ",p, pos);
            for d in 0..N_CHANNELS {
                let step_v = dir[d] * v * self.v_scale[d];
                let step_pos = pos[d] * self.step_scale[d];
                // Acceleration for correct position
                let a_pos =
                    step_pos - ((self.channels[d].pos - s0[d])
                                + 2* self.channels[d].v as i64) as f64;
                // Acceleration for correct speed
                let a_v = step_v - self.channels[d].v as f64;
                println!("a_pos: {} a_v: {}", a_pos, a_v);
                acc[d] = ((a_v + a_pos) / 2.0).round() as i16;
            }
            self.add_acc_interval(acc[X_INDEX], acc[Y_INDEX], 1);
            println!("channel: {:?} ", self.channels);
        }
    }
    
    /* Bezier curve from current position to (p2x, p2y).  Control
    points are relative to their closest endpoint. Traverses the curve
    at the given speed. */
    
    pub fn curve_to(&mut self, c1x:f64,c1y:f64, c2x:f64, c2y:f64, 
                    p2x:f64,p2y:f64, v: f64)
        {
            /*
        let x0 = self.channels[X_INDEX].pos as f64;
        let y0 = self.channels[Y_INDEX].pos as f64;
        
        let dist = (c1x*c1x + c1y*c1y).sqrt();
        let scale = if dist < 0.0001 {0.0} else {v / dist};
        let vx = c1x * scale;
        let vy = c1y * scale;

        /*println!("In: ({},{}) Out: ({},{})", 
                self.channels[X_INDEX].v, self.channels[Y_INDEX].v,
                vx, vy);*/

        let v0x = self.channels[X_INDEX].v as f64;
        let v0y = self.channels[Y_INDEX].v as f64;
        // Adjust velocities if difference is too big
        if !velocity_matches(v0x, v0y,
                             vx, vy,
                             ANGLE_SIN_LIMIT,
                             SPEED_REL_LIMIT) {
            self.goto_speed(x0, y0, vx, vy);
        }
        
        self.events.push(StepperEvent {
            ticks: u32::try_from(self.channels[X_INDEX].ticks-self.event_ticks).unwrap(),
            cmd: Command::Weight(self.curve_weight)
        });
        self.event_ticks = self.channels[X_INDEX].ticks;
        
        let x0 = self.channels[X_INDEX].pos as f64;
        let y0 = self.channels[Y_INDEX].pos as f64;
        
        let px = [x0 as f64, 
                  (x0 + c1x) as f64, 
                  (p2x + c2x) as f64, p2x as f64];
        let py = [y0 as f64, 
                  (y0 + c1y) as f64, 
                  (p2y + c2y) as f64, p2y as f64];
        let (v_adj, uvec) = bezier::find_equidistant_steps(&px,&py, v*2.0, 0.1);
        let mut acc = [Vec::<AccSegment>::new(), Vec::<AccSegment>::new()];
        let mut dvx = self.channels[X_INDEX].v; 
        let mut dvy = self.channels[Y_INDEX].v;
        /*
        let (x,y) = {
            let mut cb = |ax,ay|
            {
                acc[X_INDEX].push(AccSegment {interval: 1u16, acc: ax});
                acc[Y_INDEX].push(AccSegment {interval: 1u16, acc: ay});
                
                dvx += ax as i32;
                dvy += ay as i32;
                
            };
            bezier::constant_speed(&px,&py, v_adj*0.5, &uvec,
                                   self.channels[X_INDEX].v,
                                   self.channels[Y_INDEX].v,
                                   self.v_max, self.a_max as i16,
                                   &mut cb)
        };
        assert_eq!(y - self.channels[Y_INDEX].pos, 
                   acc[Y_INDEX].acc_distance(self.channels[Y_INDEX].v).0);
        {
            let (dx2, dv2) =
                acc[X_INDEX].acc_distance(self.channels[X_INDEX].v); 
            assert_eq!(x - self.channels[X_INDEX].pos, dx2);
            assert_eq!(dvx, dv2);

            let (dy2, dv2) =
                acc[Y_INDEX].acc_distance(self.channels[Y_INDEX].v); 
            assert_eq!(y - self.channels[Y_INDEX].pos, dy2);
            assert_eq!(dvy, dv2);

            
        }
                       
           
        
    /*println!("D: ({}, {})", x - self.channels[X_INDEX].pos,
                 y - self.channels[Y_INDEX].pos);             */
        self.channels[X_INDEX].pos = x;
        self.channels[Y_INDEX].pos = y;

        self.acc_segs_to_events(&acc);
        assert_eq!(self.channels[X_INDEX].ticks, self.channels[Y_INDEX].ticks);
    /*println!("Dist: ({}, {})", 
                 acc[X_INDEX].acc_distance(self.channels[X_INDEX].v).0,
    acc[Y_INDEX].acc_distance(self.channels[Y_INDEX].v).0);*/
        self.channels[X_INDEX].v = dvx;
        self.channels[Y_INDEX].v = dvy;

        /*
        let (s,v) = self.event_distance(self.channels[X_INDEX].ticks);
        println!("s={:?}, v={:?}", s,v);
         */
         */
             */
    }
    
    fn rotate(x:f64, y:f64, sin:f64, cos:f64) -> (f64,f64)
    {
        (x*cos - y *sin, x*sin + y*cos)
    }
    pub fn arc_to(&mut self, rx:f64,ry:f64, start:f64, end:f64, rot:f64, v: f32)
        {
            /*
        const ERR_WEIGHT:f64 = 0.2;
        let x0 = self.channels[X_INDEX].pos;
        let y0 = self.channels[Y_INDEX].pos;
        let (rs, rc) = rot.sin_cos();

        let (mut vx,mut vy) = Self::rotate(-(rx as f64) * start.sin(),
                                   (ry as f64) * start.cos(),
                                   rs,rc);
        if start > end {
            vx = -vx;
            vy = -vy;
        }
        let (vx, vy) = vel_dir(vx,vy, v);

        // Adjust velocities if difference is too big
        if !velocity_matches(self.channels[X_INDEX].v as f64,
                             self.channels[Y_INDEX].v as f64,
                             vx as f64, vy as f64,
                             ANGLE_SIN_LIMIT,
                             SPEED_REL_LIMIT) {
            self.goto_speed(x0, y0, vx, vy);
        }

        self.events.push(StepperEvent {
            ticks: (self.channels[X_INDEX].ticks-self.event_ticks) as u32,
            cmd: Command::Weight(self.curve_weight)
        });
        self.event_ticks = self.channels[X_INDEX].ticks;
        
        let x0 = self.channels[X_INDEX].pos;
        let y0 = self.channels[Y_INDEX].pos;
        
        let mut vx = self.channels[X_INDEX].v; 
        let mut vy = self.channels[Y_INDEX].v;
        
        if rx == ry {
            let rx = rx as f64;
            // It's a circle just add angles
            let start = start + rot;
            let end = end + rot;
            let arc_length = (end - start).abs() * (rx as f64);
            // Find a whole number of steps
            let n_steps = (arc_length / (2.0*v as f64)).ceil() as i32;
            let step_angle = (end - start) / n_steps as f64;

            // Start position relative circle center
            let mut x = (rx * start.cos()).round() as i64;
            let mut y = (rx * start.sin()).round() as i64;

            // Move to center
            let x0 = x0 - x;
            let y0 = y0 - y;
            
            let radial_acc = rx * step_angle * step_angle*0.5;
            
            let mut acc = [Vec::<AccSegment>::new(), Vec::<AccSegment>::new()];
            for step in 0 .. n_steps {
                let angle = step as f64 * step_angle + start;
                let x1 = (rx * (angle + step_angle).cos()).round() as i64;
                let y1 = (rx * (angle + step_angle).sin()).round() as i64;

                let ax = -radial_acc * (angle+ step_angle/2.0).cos();
                let ax = (ax + ((x1 - x) as f64 - (2*vx) as f64 - ax) 
                          * ERR_WEIGHT).round() as i16;
                
                let ay = -radial_acc* (angle+step_angle/2.0).sin();
                let ay = (ay + ((y1 - y) as f64 - (2*vy) as f64 - ay) 
                          * ERR_WEIGHT).round() as i16;

                acc[X_INDEX].push(AccSegment {interval: 1u16,
                                              acc: ax});
                acc[Y_INDEX].push(AccSegment {interval: 1u16,
                                              acc: ay});
                
                x += 2*vx as i64 + ax as i64;
                y += 2*vy as i64 + ay as i64;
                vx += ax as i32;
                vy += ay as i32;
            }
            self.acc_segs_to_events(&acc);
            self.channels[X_INDEX].pos = x + x0;
            self.channels[Y_INDEX].pos = y + y0;

            self.channels[X_INDEX].v = vx;
            self.channels[Y_INDEX].v = vy;

        } else {
            panic!("Non circular ellipses not supported");
            // TODO
        }
             */
    }
    
    pub fn set_weight(&mut self, weight: i32)
    {
        self.curve_weight = weight;
    }


    
    pub fn draw_curves(&mut self, segs: &[CurveSegment], v: f64)
        {
            /*
        let mut iter = segs.iter();
        let mut seg =
            match iter.next() {
                Some(s) => s,
                None => return
            };
        loop {
            // Start points
            let x0 = self.channels[X_INDEX].pos;
            let y0 = self.channels[Y_INDEX].pos;
            let (end_x, end_y) = match *seg {
                CurveSegment::GoTo(x,y) => (x,y), 
                CurveSegment::GoToRel(x,y) =>  (x + x0, y + y0),
                CurveSegment::LineTo(x,y) => (x,y),
                CurveSegment::LineToRel(x,y) =>  (x + x0, y + y0),
                CurveSegment::CurveTo(x, y, _, _, _, _) => (x, y), 
                CurveSegment::CurveToRel(x, y, _, _, _, _) => (x + x0, y + y0),
                CurveSegment::Arc(rx,ry, start,end, rot) => {
                    let x = (rx as f64)*(end.cos()-start.cos());
                    let y = (ry as f64)*(end.sin()-start.sin());
                    let s = rot.sin();
                    let c = rot.cos();
                    ((x*c - y *s) as i64 + x0,  (x*s + y*c) as i64 + y0)
                }
            };
            
            let next_seg = iter.next();
            //println!("({},{}) -> ({},{})", x0,y0, end_x, end_y);
            // Skip zero length gotos and lines
            if match *seg {
                // Curves may start and end in the same point as they start
                CurveSegment::CurveTo(_,_,_,_,_,_) => true,
                CurveSegment::CurveToRel(_,_,_,_,_,_) => true,
                _ => {
                    // Skip lines and gotos that are too short
                    let x = end_x - x0;
                    let y = end_y - y0;
                    (x*x + y*y) >= (self.a_max / 10) as i64
                }
            }
            {
                let next_speed: Option<(i32,i32)>;
                if let Some(next_seg) = next_seg {
                    // Calculate initial velocity for next segment
                    next_speed = match *next_seg {
                        CurveSegment::GoTo(_,_) 
                            | CurveSegment::GoToRel(_,_) => None,
                        CurveSegment::LineTo(x,y) => {
                            if (x - end_x) == 0 && (y - end_y) == 0 {
                                None
                            } else {
                                Some(vel_dir((x - end_x) as f64,
                                             (y - end_y) as f64, v))
                            }
                        },
                        CurveSegment::LineToRel(x,y) => {
                            if x == 0 && y == 0 {
                                None
                            } else {
                                Some(vel_dir(x as f64,y as f64,v))
                            }
                        },
                        CurveSegment::CurveTo(_, _, c1x, c1y, _, _) 
                            | CurveSegment::CurveToRel(_, _, c1x, c1y, _, _)=> {
                                if c1x == 0 && c1y == 0 {
                                    None
                                } else {
                                    Some(vel_dir(c1x as f64,c1y as f64,v))
                                }
                            },
                        CurveSegment::Arc(rx,ry, start,end, rot) => {
                            let mut vx = -(rx as f64) * start.sin();
                            let mut vy = (ry as f64) * start.cos();
                            if start > end {
                                vx = -vx;
                                vy = -vy;
                            }
                            let (s,c) = rot.sin_cos();
                            Some(vel_dir(vx*c - vy *s,
                                         vx*s + vy*c, v))
                        }
                    };
                } else {
                    next_speed = Some((0,0));
                }
                println!("Seg: {:?}, next {:?} {:?}", seg, next_seg, next_speed);
                match (seg, next_speed) {
                    (&CurveSegment::GoTo(x,y), None) =>
                        self.goto_speed(x,y,0,0),
                    (&CurveSegment::GoTo(x,y), Some((vx,vy))) => 
                        self.goto_speed(x,y,vx,vy),
                    
                    (&CurveSegment::GoToRel(x,y), None) => {
                        self.goto_speed(x + x0, y + y0,0,0);
                    },
                    (&CurveSegment::GoToRel(x,y), Some((vx,vy))) => {
                        self.goto_speed(x + x0, y + y0,vx,vy);
                    },
                    
                    (&CurveSegment::LineTo(x,y), _) => {
                        let cx = (x - x0) / 3;
                        let cy = (y - y0) / 3;
                        self.curve_to(cx, cy, -cx, -cy, x, y, v);
                    },
                    
                    (&CurveSegment::LineToRel(x, y), _) => {
                        let cx = x / 3;
                        let cy = y / 3;
                        self.curve_to(cx, cy, -cx, -cy, x + x0, y + y0, v);
                    },

                    
                    (&CurveSegment::CurveTo(x,y,c1x, c1y, c2x, c2y), _) => {
                        self.curve_to(c1x, c1y, c2x, c2y, x, y, v);
                    },
                    
                    (&CurveSegment::CurveToRel(x,y,c1x,c1y, c2x,c2y), _) => {
                        self.curve_to(c1x, c1y, c2x, c2y, x + x0, y + y0, v);
                        
                    },
                    (&CurveSegment::Arc(rx,ry,start,end, rot), _) => {
                        self.arc_to(rx,ry,start,end, rot, v);
                    }
                }
            }
            if let Some(s) = next_seg {
                seg = s;
            } else {
                break;
            }
        }
        */
        
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

    let line = Line::new(23.0,45.0);
    ctxt.approx_curve(&line, 1.0);

}

#[test]
fn test_approx_curve_ellipse()
{
    let mut ctxt = StepperContext::new(&[10, 10],
                                       &[10, 10],
                                       &[3.0, 7.0],
                                       &[1.5, 3.5]);

    let circle = Ellipse::new(670.0,670.0, 0.0, PI/2.0, 0.0);
    ctxt.approx_curve(&circle, 30.0);

}
