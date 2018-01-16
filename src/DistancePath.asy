void draw_axis(real taxismax, 
	       real vaxismin, 	       
	       real vaxismax) {
  draw((0,vaxismin) -- (0,vaxismax), arrow=Arrow);
  draw((0,0) -- (taxismax,0), arrow=Arrow);
}

void draw_dist(pair a, pair b, Label value, pair offset) {
  real c = dot((b-a),offset) / dot(offset , offset);
  pair bo;
  pair ao;
  if (c >= 0) {
    bo = b+offset;
    ao = a+offset * (1 + c);
  } else {
    bo = b+offset * (1 - c);
    ao = a+offset;
  }
  draw(b--bo ^^ a--ao);
  Label L = Label(value, position = MidPoint);
  draw(ao -- bo, arrow=Arrows, L=L);
  
}
  
size(10cm);

real taxismax = 15;
real vaxismax = 10;
real vaxismin = -1;


draw_axis(taxismax, vaxismin, vaxismax);

real a = 3;
real v0 = 1;
real vn = 3;
real vmax = 9;
real dt1 = 2;
real dt2 = 9;
real dt3 = 3;
real dt4 = 1;

real v1 = v0+a*dt1;
real t1 = dt1;

real v2 = vmax;
real t2 = t1 + 1;

real v3 = vmax;
real t3 = dt3 + t1;

real v4 = vmax - 1;
real t4 = t3 + 1;

real v5 = v4;
real t5 = t1 + dt2-1;

real v6 = a*dt4 + vn;
real t6 = t1 + dt2;

real v7 = vn;
real t7 = t6 + dt4;


draw((0,v0) -- (t1,v1) -- (t2,v2) -- (t3,v3) -- (t4, v4) -- (t5, v5) -- (t6,v6) -- (t7,v7), red+linewidth(1pt));
label(Label("$\mathtt{v0}$", align= LeftSide), position=(0,v0), black);
label(Label("$\mathtt{vn}$", align= LeftSide), position=(0,vn), black);
label(Label("$\mathtt{vflat}$", align= LeftSide), position=(0,vmax), black);

draw_dist((0,v0),(t7,vn), Label("$\mathtt{t\_total}$",align=0.5N), (0,-0.5)); 
draw_dist((t3,v3),(t5,v5), Label("$|\mathtt{t\_adjust}|$",align=0.5N), (0,0.5)); 

