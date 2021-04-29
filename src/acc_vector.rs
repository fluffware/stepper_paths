use std::slice;
use std::fmt;

/// Acceleration for a given number of steps
pub struct AccSegment
{
    /// Number of steps to accelerate for
    pub interval:u16,
    /// Acceleration
    pub acc:i16
}

impl fmt::Debug for AccSegment {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        if self.interval == 1 {
            write!(f, "{:?}", self.acc)
        } else {
            write!(f, "{:?} *{:?}", self.acc, self.interval)
        }
    }
}

pub trait AccVector<'a> {
    type I;
    fn acc_iter(&'a self, t: u16) -> Self::I;
    fn acc_push(&'a mut self, t: u16, a: i16);
    fn acc_distance(&'a self, v0: i32) -> (i64,i32); // Returns (dist, velocity)
}

pub struct AccIter<'a, I> 
    where I : Iterator<Item = &'a AccSegment> + Sized {
    count: u16,
    acc: i16,
    seg: I
}

impl<'a, I> Iterator for AccIter<'a, I>
    where I : Iterator<Item = &'a AccSegment>
{
    type Item = i16;
    fn next(&mut self) -> Option<i16> {
        self.count -= 1;
        if self.count == 0 {
            match self.seg.next() {
                Some(&AccSegment{interval, acc}) => {
                    self.count = interval; 
                    self.acc = acc;
                },
                None => return None
            }
        }
        Some(self.acc)
        
    }
}


impl<'a> AccVector<'a> for Vec<AccSegment>
{
    type I= AccIter<'a, slice::Iter<'a, AccSegment>>;
    fn acc_iter(&'a self, t: u16) -> Self::I
    {
        let mut start: u16 = 0;
        let mut seg = self.iter();
        while let Some(&AccSegment{interval: i,acc: a}) = seg.next() {
            let end = start + i;
            if t >= start && t < end {
                return AccIter{count:end - t + 1, acc:a , seg};
            }
            start = end;
        }
        panic!("Time outside range");
    }

    fn acc_push(&'a mut self, t: u16, a: i16)
    {
        if t == 0 {return;}
        if match self.last_mut() {
            Some(&mut AccSegment{ref mut interval, acc}) => {
                if acc == a {
                    *interval += t;
                    false                    
                } else {
                    true
                }
                
            },
            None => true
        } {
            self.push(AccSegment{interval: t, acc: a});
        }
    }
    
    fn acc_distance(&'a self, v0:i32) -> (i64, i32) {
        let mut s = 0i64;
        let mut v = v0 as i32;
        for &AccSegment{interval: i,acc: a} in self.iter() {
            let v_next = v + a as i32 * i as i32;
            s += (v + v_next) as i64 * i as i64;
            v = v_next;
        }
        (s,v)
    }
}


#[test]
fn test_acc_iter()
{
    let v = vec![AccSegment{interval: 2, acc:8},
                 AccSegment{interval: 1, acc:4},
                 AccSegment{interval: 3, acc:-7}];
    let i = AccIter{count:1, acc:9, seg: v.iter()};
    let v2: Vec<i16> = i.collect();
    assert_eq!(v2, vec![8,8,4,-7,-7,-7]);
}

#[test]
fn test_acc_vector()
{
    let v = vec![AccSegment{interval: 2, acc:8},
                 AccSegment{interval: 1, acc:4},
                 AccSegment{interval: 3, acc:-7}];
    assert_eq!(v.acc_iter(0).collect::<Vec<i16>>(), vec![8,8,4,-7,-7,-7]);
    assert_eq!(v.acc_iter(2).collect::<Vec<i16>>(), vec![4,-7,-7,-7]);
    assert_eq!(v.acc_iter(4).collect::<Vec<i16>>(), vec![-7,-7]);
    assert_eq!(v.acc_iter(5).collect::<Vec<i16>>(), vec![-7]);

   
}
