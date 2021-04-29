extern crate libc;
extern crate core;
extern crate serial;
use stepper_context::StepperEvent;
use stepper_context::Command;
use std;

use std::io::Write;
use std::io::Read;
use std::io::ErrorKind;
use std::io;
use std::thread;
use std::time;

fn read_reply(inp: &mut dyn Read) -> io::Result<String> {
    let mut reply = String::new();
    let mut buf = [0u8;1];
    loop {
        match inp.read(& mut buf) {
            Ok(s) => {
                if s == 1 {
                    if buf[0] == b'\n' {
                        //println!("Reply: '{}'", reply);
                        return Ok(reply);
                    } else {
                        reply.push(buf[0] as char);
                    }
                }
            }
            Err(ref e) if e.kind() == std::io::ErrorKind::TimedOut => {},
            Err(e) => return Err(e)
        };
    }

}

fn request_time<T: Read+Write>(serport: &mut T) -> io::Result<i64> {
    serport.write_all(b"T\n").unwrap();
    read_reply(serport)
        .and_then(|r| i64::from_str_radix(&r,10)
                  .map_err(|_| {
                      io::Error::new(ErrorKind::Other,
                                     "Invalid integer in reply")
                  }))
}


fn send_cmds<T>(serport: &mut T, t0:i64, events: &[StepperEvent])
    where T: Read + Write
{
    let mut t = t0;
    let mut cmds = String::new();
    for ev in events {
        t += ev.ticks as i64;
        match ev.cmd {
            Command::Acc(ch, acc) => {
                cmds += &format!("A{},{},{}\n", ch, t, acc);
            },
            Command::Weight(v) => {
                cmds += &format!("L{},{}\n", t, v);
            }
        }
    }
    
    serport.write_all(cmds.as_bytes()).unwrap();
    //println!("Cmds: {}", cmds);
}

const TICKS_PER_SEC:i64 = 128;
pub fn play_events<T>(serport: &mut T, events: &[StepperEvent]) -> io::Result<()>
    where T: Read + Write
{
    let mut t = request_time(serport).unwrap() + TICKS_PER_SEC;
    let mut pending = Vec::<StepperEvent>::new();
    let mut ev_iter = events.iter();
    loop {
        while pending.is_empty() {
            match ev_iter.next() {
                Some(e) =>
                    pending.push(e.clone()),
                None => {break;}
            }
        }
        if pending.is_empty() {break;}
        send_cmds(serport, t, &pending);
        let mut first_full : i32 = -1;
        for (i,p) in pending.iter().enumerate() {
            match read_reply(serport) {
                Ok(r) => {
                    //println!("'{}'", r);
                    if r != "OK" {
                        if r == "ERR 1" { // Queue is full
                            //println!("Queue full");
                            if first_full < 0 {
                                first_full = i as i32;
                            }
                        } else {
                            serport.write_all("R\n".as_bytes()).unwrap();    
                            return Err(
                                io::Error::new(ErrorKind::Other, 
                                               format!("Contoller error: {}",r)));
                        }
                    }
                    if first_full == -1 {
                        t += p.ticks as i64;
                    }
                },
                Err(e) => return  Err(e)
            }
        }
        if first_full == -1 {
            pending.clear();
        } else {
            pending.drain(0..(first_full as usize));
            thread::sleep(time::Duration::from_millis(10));
        }
    }
    Ok(())

}

