#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {perspective
  right -10.55*x up 10.55*y
  direction 50.00*z
  location <0,0,50.00> look_at <0,0,0>}


light_source {<  2.00,   3.00,  40.00> color White
  area_light <0.70, 0, 0>, <0, 0.70, 0>, 3, 3
  adaptive 1 jitter}
// no fog
#declare simple = finish {phong 0.7}
#declare pale = finish {ambient 0.5 diffuse 0.85 roughness 0.001 specular 0.200 }
#declare intermediate = finish {ambient 0.3 diffuse 0.6 specular 0.1 roughness 0.04}
#declare vmd = finish {ambient 0.0 diffuse 0.65 phong 0.1 phong_size 40.0 specular 0.5 }
#declare jmol = finish {ambient 0.2 diffuse 0.6 specular 1 roughness 0.001 metallic}
#declare ase2 = finish {ambient 0.05 brilliance 3 diffuse 0.6 metallic specular 0.7 roughness 0.04 reflection 0.15}
#declare ase3 = finish {ambient 0.15 brilliance 2 diffuse 0.6 metallic specular 1.0 roughness 0.001 reflection 0.0}
#declare glass = finish {ambient 0.05 diffuse 0.3 specular 1.0 roughness 0.001}
#declare glass2 = finish {ambient 0.01 diffuse 0.3 specular 1.0 reflection 0.25 roughness 0.001}
#declare Rcell = 0.050;
#declare Rbond = 0.070;

#macro atom(LOC, R, COL, TRANS, FIN)
  sphere{LOC, R texture{pigment{color COL transmit TRANS} finish{FIN}}}
#end
#macro constrain(LOC, R, COL, TRANS FIN)
union{torus{R, Rcell rotate 45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     torus{R, Rcell rotate -45*z texture{pigment{color COL transmit TRANS} finish{FIN}}}
     translate LOC}
#end

cylinder {< -4.52,  -4.52,  -0.00>, < -4.52,  -4.52,  -9.04>, Rcell pigment {Black}}
cylinder {< -4.52,   4.52,  -0.00>, < -4.52,   4.52,  -9.04>, Rcell pigment {Black}}
cylinder {<  4.52,   4.52,   0.00>, <  4.52,   4.52,  -9.04>, Rcell pigment {Black}}
cylinder {<  4.52,  -4.52,   0.00>, <  4.52,  -4.52,  -9.04>, Rcell pigment {Black}}
cylinder {< -4.52,  -4.52,  -0.00>, < -4.52,   4.52,  -0.00>, Rcell pigment {Black}}
cylinder {< -4.52,  -4.52,  -9.04>, < -4.52,   4.52,  -9.04>, Rcell pigment {Black}}
cylinder {<  4.52,  -4.52,  -9.04>, <  4.52,   4.52,  -9.04>, Rcell pigment {Black}}
cylinder {<  4.52,  -4.52,   0.00>, <  4.52,   4.52,   0.00>, Rcell pigment {Black}}
cylinder {< -4.52,  -4.52,  -0.00>, <  4.52,  -4.52,   0.00>, Rcell pigment {Black}}
cylinder {< -4.52,  -4.52,  -9.04>, <  4.52,  -4.52,  -9.04>, Rcell pigment {Black}}
cylinder {< -4.52,   4.52,  -9.04>, <  4.52,   4.52,  -9.04>, Rcell pigment {Black}}
cylinder {< -4.52,   4.52,  -0.00>, <  4.52,   4.52,   0.00>, Rcell pigment {Black}}
atom(< -4.52,  -4.52,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #0
atom(< -4.52,   0.00,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #1
atom(<  0.00,  -4.52,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #2
atom(<  0.00,   0.00,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #3
atom(<  0.00,   0.00,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #4
atom(<  0.00,  -4.52,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #5
atom(< -4.52,   0.00,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #6
atom(< -4.52,  -4.52,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #7
atom(<  3.69,  -3.75,  -6.78>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #8
atom(<  3.69,  -0.77,  -6.78>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #9
atom(<  0.83,  -3.75,  -6.78>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #10
atom(<  0.83,  -0.77,  -6.78>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #11
atom(< -3.75,   2.26,  -8.21>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #12
atom(< -0.77,   2.26,  -8.21>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #13
atom(< -3.75,   2.26,  -5.35>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #14
atom(< -0.77,   2.26,  -5.35>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #15
atom(<  2.26,   3.69,  -0.77>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #16
atom(<  2.26,   3.69,  -3.75>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #17
atom(<  2.26,   0.83,  -0.77>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #18
atom(<  2.26,   0.83,  -3.75>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #19
atom(< -0.83,   0.77,  -2.26>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #20
atom(< -0.83,   3.75,  -2.26>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #21
atom(< -3.69,   0.77,  -2.26>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #22
atom(< -3.69,   3.75,  -2.26>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #23
atom(<  0.77,  -2.26,  -3.69>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #24
atom(<  3.75,  -2.26,  -3.69>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #25
atom(<  0.77,  -2.26,  -0.83>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #26
atom(<  3.75,  -2.26,  -0.83>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #27
atom(< -2.26,  -0.83,  -5.29>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #28
atom(< -2.26,  -0.83,  -8.27>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #29
atom(< -2.26,  -3.69,  -5.29>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #30
atom(< -2.26,  -3.69,  -8.27>, 0.56, rgb <0.62, 0.39, 0.71>, 0.0, ase3) // #31
atom(<  4.52,  -4.52,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #32
atom(<  4.52,   0.00,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #33
atom(<  4.52,   0.00,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #34
atom(<  4.52,  -4.52,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #35
atom(< -4.52,   4.52,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #36
atom(<  0.00,   4.52,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #37
atom(<  0.00,   4.52,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #38
atom(< -4.52,   4.52,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #39
atom(<  4.52,   4.52,  -0.00>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #40
atom(<  4.52,   4.52,  -4.52>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #41
atom(< -4.52,  -4.52,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #42
atom(<  0.00,   0.00,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #43
atom(<  0.00,  -4.52,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #44
atom(< -4.52,   0.00,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #45
atom(<  4.52,  -4.52,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #46
atom(<  4.52,   0.00,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #47
atom(< -4.52,   4.52,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #48
atom(<  0.00,   4.52,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #49
atom(<  4.52,   4.52,  -9.04>, 0.50, rgb <0.94, 0.56, 0.63>, 0.0, ase3) // #50
cylinder {<  0.00,   0.00,  -4.52>, < -2.26,   0.00,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -4.52>, < -2.26,   0.00,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -4.52>, <  0.00,  -2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -4.52>, <  0.00,  -2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -4.52>, <  0.00,   0.00,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -0.00>, <  0.00,   0.00,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -0.00>, < -2.26,  -4.52,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -0.00>, < -2.26,  -4.52,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -0.00>, <  0.00,  -4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -4.52>, <  0.00,  -4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -0.00>, <  0.00,  -2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -0.00>, <  0.00,  -2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -0.00>, < -4.52,  -2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -0.00>, < -4.52,  -2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -0.00>, < -4.52,   0.00,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -4.52>, < -4.52,   0.00,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -0.00>, < -2.26,   0.00,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -0.00>, < -2.26,   0.00,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -4.52>, < -4.52,  -4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -0.00>, < -4.52,  -4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -4.52>, < -4.52,  -2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -4.52>, < -4.52,  -2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -4.52>, < -2.26,  -4.52,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -4.52>, < -2.26,  -4.52,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  3.69,  -0.77,  -6.78>, <  3.69,  -2.26,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  3.69,  -3.75,  -6.78>, <  3.69,  -2.26,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  0.83,  -3.75,  -6.78>, <  2.26,  -3.75,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  3.69,  -3.75,  -6.78>, <  2.26,  -3.75,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  0.83,  -0.77,  -6.78>, <  2.26,  -0.77,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  3.69,  -0.77,  -6.78>, <  2.26,  -0.77,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  0.83,  -0.77,  -6.78>, <  0.83,  -2.26,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  0.83,  -3.75,  -6.78>, <  0.83,  -2.26,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -0.77,   2.26,  -8.21>, < -2.26,   2.26,  -8.21>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -3.75,   2.26,  -8.21>, < -2.26,   2.26,  -8.21>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -3.75,   2.26,  -5.35>, < -3.75,   2.26,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -3.75,   2.26,  -8.21>, < -3.75,   2.26,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -0.77,   2.26,  -5.35>, < -0.77,   2.26,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -0.77,   2.26,  -8.21>, < -0.77,   2.26,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -0.77,   2.26,  -5.35>, < -2.26,   2.26,  -5.35>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -3.75,   2.26,  -5.35>, < -2.26,   2.26,  -5.35>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  2.26,   3.69,  -3.75>, <  2.26,   3.69,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  2.26,   3.69,  -0.77>, <  2.26,   3.69,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  2.26,   0.83,  -0.77>, <  2.26,   2.26,  -0.77>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  2.26,   3.69,  -0.77>, <  2.26,   2.26,  -0.77>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  2.26,   0.83,  -3.75>, <  2.26,   2.26,  -3.75>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  2.26,   3.69,  -3.75>, <  2.26,   2.26,  -3.75>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  2.26,   0.83,  -3.75>, <  2.26,   0.83,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  2.26,   0.83,  -0.77>, <  2.26,   0.83,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -0.83,   3.75,  -2.26>, < -0.83,   2.26,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -0.83,   0.77,  -2.26>, < -0.83,   2.26,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -3.69,   0.77,  -2.26>, < -2.26,   0.77,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -0.83,   0.77,  -2.26>, < -2.26,   0.77,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -3.69,   3.75,  -2.26>, < -2.26,   3.75,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -0.83,   3.75,  -2.26>, < -2.26,   3.75,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -3.69,   3.75,  -2.26>, < -3.69,   2.26,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -3.69,   0.77,  -2.26>, < -3.69,   2.26,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  3.75,  -2.26,  -3.69>, <  2.26,  -2.26,  -3.69>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  0.77,  -2.26,  -3.69>, <  2.26,  -2.26,  -3.69>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  0.77,  -2.26,  -0.83>, <  0.77,  -2.26,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  0.77,  -2.26,  -3.69>, <  0.77,  -2.26,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  3.75,  -2.26,  -0.83>, <  3.75,  -2.26,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  3.75,  -2.26,  -3.69>, <  3.75,  -2.26,  -2.26>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  3.75,  -2.26,  -0.83>, <  2.26,  -2.26,  -0.83>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  0.77,  -2.26,  -0.83>, <  2.26,  -2.26,  -0.83>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -2.26,  -0.83,  -8.27>, < -2.26,  -0.83,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -2.26,  -0.83,  -5.29>, < -2.26,  -0.83,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -2.26,  -3.69,  -5.29>, < -2.26,  -2.26,  -5.29>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -2.26,  -0.83,  -5.29>, < -2.26,  -2.26,  -5.29>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -2.26,  -3.69,  -8.27>, < -2.26,  -2.26,  -8.27>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -2.26,  -0.83,  -8.27>, < -2.26,  -2.26,  -8.27>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -2.26,  -3.69,  -8.27>, < -2.26,  -3.69,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {< -2.26,  -3.69,  -5.29>, < -2.26,  -3.69,  -6.78>, Rbond texture{pigment {color rgb <0.62, 0.39, 0.71> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -0.00>, <  2.26,  -4.52,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -0.00>, <  2.26,  -4.52,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -4.52>, <  2.26,   0.00,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -4.52>, <  2.26,   0.00,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -0.00>, <  2.26,   0.00,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -0.00>, <  2.26,   0.00,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -0.00>, <  4.52,  -2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -0.00>, <  4.52,  -2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -0.00>, <  4.52,   0.00,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -4.52>, <  4.52,   0.00,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -4.52>, <  2.26,  -4.52,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -4.52>, <  2.26,  -4.52,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -4.52>, <  4.52,  -4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -0.00>, <  4.52,  -4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -4.52>, <  4.52,  -2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -4.52>, <  4.52,  -2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -0.00>, < -4.52,   2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -0.00>, < -4.52,   2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -4.52>, <  0.00,   2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -4.52>, <  0.00,   2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -0.00>, <  0.00,   2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -0.00>, <  0.00,   2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -0.00>, < -2.26,   4.52,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -0.00>, < -2.26,   4.52,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -0.00>, <  0.00,   4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -4.52>, <  0.00,   4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -4.52>, < -4.52,   2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -4.52>, < -4.52,   2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -4.52>, < -4.52,   4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -0.00>, < -4.52,   4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -4.52>, < -2.26,   4.52,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -4.52>, < -2.26,   4.52,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -0.00>, <  4.52,   2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -0.00>, <  4.52,   2.26,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -0.00>, <  2.26,   4.52,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -0.00>, <  2.26,   4.52,  -0.00>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -4.52>, <  4.52,   2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -4.52>, <  4.52,   2.26,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -4.52>, <  2.26,   4.52,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -4.52>, <  2.26,   4.52,  -4.52>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -4.52>, <  4.52,   4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -0.00>, <  4.52,   4.52,  -2.26>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -9.04>, < -4.52,  -4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -4.52>, < -4.52,  -4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -9.04>, <  0.00,   0.00,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -4.52>, <  0.00,   0.00,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -9.04>, <  0.00,  -4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -4.52>, <  0.00,  -4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -9.04>, < -2.26,  -4.52,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -9.04>, < -2.26,  -4.52,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -9.04>, <  0.00,  -2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -9.04>, <  0.00,  -2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -9.04>, < -4.52,   0.00,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -4.52>, < -4.52,   0.00,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -9.04>, < -4.52,  -2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,  -4.52,  -9.04>, < -4.52,  -2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -9.04>, < -2.26,   0.00,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -9.04>, < -2.26,   0.00,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -9.04>, <  4.52,  -4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -4.52>, <  4.52,  -4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -9.04>, <  2.26,  -4.52,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -4.52,  -9.04>, <  2.26,  -4.52,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -9.04>, <  4.52,   0.00,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -4.52>, <  4.52,   0.00,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -9.04>, <  2.26,   0.00,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -9.04>, <  2.26,   0.00,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -9.04>, <  4.52,  -2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,  -4.52,  -9.04>, <  4.52,  -2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -9.04>, < -4.52,   4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -4.52>, < -4.52,   4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -9.04>, < -4.52,   2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   0.00,  -9.04>, < -4.52,   2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -9.04>, <  0.00,   4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -4.52>, <  0.00,   4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -9.04>, <  0.00,   2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   0.00,  -9.04>, <  0.00,   2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -9.04>, < -2.26,   4.52,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {< -4.52,   4.52,  -9.04>, < -2.26,   4.52,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -9.04>, <  4.52,   4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -4.52>, <  4.52,   4.52,  -6.78>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -9.04>, <  4.52,   2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   0.00,  -9.04>, <  4.52,   2.26,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  4.52,   4.52,  -9.04>, <  2.26,   4.52,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   4.52,  -9.04>, <  2.26,   4.52,  -9.04>, Rbond texture{pigment {color rgb <0.94, 0.56, 0.63> transmit 0.0} finish{ase3}}}
// no constraints
