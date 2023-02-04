use std::i64;

fn same_sign(a: i64, b: i64) -> bool {
    if a >= 0 {
        b >= 0
    } else {
        b < 0
    }
}

pub fn sqpoly(a: i64, b: i64, c: i64, x: i64) -> i64 {
    a * x * x + b * x - c
}

/* Find the value of x that satisfies a*x^2 + b*x >= c
Such that x+1 or x-1 will not satisfy the comparasion */

pub fn find_root(a: i64, b: i64, c: i64, min: i64, max: i64) -> Option<i64> {
    let mut minx = min;
    let mut maxx = max;
    let mut miny = sqpoly(a, b, c, minx);
    let mut maxy = sqpoly(a, b, c, maxx);
    //println!("Max: {} Min: {}", maxy, miny);
    if miny == 0 {
        return Some(minx);
    };
    if maxy == 0 {
        return Some(maxx);
    };
    if same_sign(miny, maxy) {
        return None;
    }
    loop {
        if miny == 0 {
            return Some(minx);
        };
        if maxy == 0 {
            return Some(maxx);
        };
        if maxx - minx <= 1 {
            return if maxy >= 0 { Some(maxx) } else { Some(minx) };
        }
        let midx = (minx + maxx) / 2;
        let midy = sqpoly(a, b, c, midx);
        //println!("{} => {}", midx, midy);
        if same_sign(maxy, midy) {
            maxy = midy;
            maxx = midx;
        } else {
            miny = midy;
            minx = midx;
        }
    }
}

// a > 0
pub fn find_roots(a: i64, b: i64, c: i64) -> Option<(i64, i64)> {
    let s = b * b + 4 * a * c;
    if s < 0 {
        return None;
    } // No solution
      // Find minima
    let xmin =
    // Round -b / (2*a) to neareast integer
        if b >= 0 {
            -((b + a) / (2*a))
        } else {
            (-b + a)  / (2*a)
        };

    //println!("xmin: {}", xmin);
    let ymin = sqpoly(a, b, c, xmin);
    if s == 0 && ymin == 0 {
        return Some((xmin, xmin)); // One exact solution
    }
    if ymin >= 0 {
        // Both roots are between two adjacent integers
        if -b < 2 * a * xmin {
            return Some((xmin - 1, xmin));
        } else {
            return Some((xmin, xmin + 1));
        }
    }
    let mut limit = 16;
    while sqpoly(a, b, c, xmin + limit) <= 0 {
        limit *= 2;
    }
    limit += 1;
    match (
        find_root(a, b, c, xmin - limit, xmin),
        find_root(a, b, c, xmin, xmin + limit),
    ) {
        (Some(x1), Some(x2)) => Some((x1, x2)),
        _ => None,
    }
}

#[cfg(test)]
fn ceil_idiv(n: i64, d: i64) -> i64 {
    if n == 0 {
        return 0;
    }
    if d > 0 {
        if n >= 0 {
            (n + d - 1) / d
        } else {
            -((-n) / d)
        }
    } else {
        if n >= 0 {
            -(n / -d)
        } else {
            (-n - d - 1) / -d
        }
    }
}

#[cfg(test)]
fn floor_idiv(n: i64, d: i64) -> i64 {
    if n == 0 {
        return 0;
    }
    if d > 0 {
        if n >= 0 {
            n / d
        } else {
            -((-n + d - 1) / d)
        }
    } else {
        if n >= 0 {
            -((n - d - 1) / -d)
        } else {
            (-n) / -d
        }
    }
}

#[cfg(test)]
fn round_idiv(n: i64, d: i64) -> i64 {
    if d > 0 {
        if n >= 0 {
            (n + d / 2) / d
        } else {
            -((-n + d / 2) / d)
        }
    } else {
        if n >= 0 {
            -((n + -d / 2) / -d)
        } else {
            (-n + -d / 2) / -d
        }
    }
}

#[test]
fn idiv_test() {
    for n in -13..14 {
        for d in -5..5 {
            if d == 0 {
                continue;
            }
            let ceil_q = ceil_idiv(n, d);
            let floor_q = floor_idiv(n, d);
            let round_q = round_idiv(n, d);
            println!("{}/{} = ({} <= {} <= {})", n, d, floor_q, round_q, ceil_q);
            if d > 0 {
                assert!(ceil_q * d >= n);
                assert!(floor_q * d <= n);
            } else {
                assert!(ceil_q * d <= n);
                assert!(floor_q * d >= n);
            }
            assert!((ceil_q - floor_q).abs() <= 1);
            assert!(round_q <= ceil_q);
            assert!(round_q >= floor_q);
            assert!((2 * (round_q * d - n)).abs() <= d.abs());
        }
    }
}

#[cfg(test)]
fn find_root_check_limits(a: i64, b: i64, c: i64, min: i64, max: i64) {
    let xmin = round_idiv(-b, 2 * a);
    let res = find_root(a, b, c, min, max);
    if 4 * a * c + b * b < 0 {
        assert_eq!(res, None);
        println!("No root");
    } else {
        match res {
            Some(root) => {
                println!("x = {}", root);
                assert!(sqpoly(a, b, c, root) >= 0);
                if root > xmin {
                    assert!(sqpoly(a, b, c, root - 1) < 0);
                } else if root < xmin {
                    assert!(sqpoly(a, b, c, root + 1) < 0);
                }
            }
            None => {
                assert!(sqpoly(a, b, c, xmin) >= 0);
            }
        }
    }
}

#[cfg(test)]
fn find_root_check(a: i64, b: i64, c: i64) {
    const MAX_LIMIT: i64 = 2642245;
    const MIN_LIMIT: i64 = -2642245;

    print!("{}*x^2 + {}*x = {} => ", a, b, c);
    let xmin = round_idiv(-b, 2 * a);
    find_root_check_limits(a, b, c, xmin, MAX_LIMIT);
    find_root_check_limits(a, b, c, MIN_LIMIT, xmin);
}

#[test]
fn find_root_test() {
    for a in 1..3 {
        for b in -5..5 {
            for c in -1..64 {
                find_root_check(a, b, c);
            }
        }
    }
}

#[test]
fn find_roots_test() {
    assert_eq!(find_roots(4, -12, -8), Some((1, 2)));
    assert_eq!(find_roots(4, -12, -7), Some((0, 3)));
    assert_eq!(find_roots(9, -90, -225), Some((5, 5)));
    assert_eq!(find_roots(4, -12, -9), Some((1, 2)));
    assert_eq!(find_roots(9, -27, -20), Some((1, 2)));
    assert_eq!(find_roots(9, -21, -12), Some((1, 2)));
}
