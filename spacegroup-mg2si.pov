#include "colors.inc"
#include "finish.inc"

global_settings {assumed_gamma 1 max_trace_level 6}
background {color White}
camera {perspective
  right -9.11*x up 8.93*y
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

cylinder {< -3.89,  -3.52,  -2.73>, < -2.25,  -2.46,  -8.78>, Rcell pigment {Black}}
cylinder {< -3.89,   2.74,  -1.62>, < -2.25,   3.81,  -7.67>, Rcell pigment {Black}}
cylinder {<  2.25,   2.46,   0.00>, <  3.89,   3.52,  -6.05>, Rcell pigment {Black}}
cylinder {<  2.25,  -3.81,  -1.10>, <  3.89,  -2.74,  -7.15>, Rcell pigment {Black}}
cylinder {< -3.89,  -3.52,  -2.73>, < -3.89,   2.74,  -1.62>, Rcell pigment {Black}}
cylinder {< -2.25,  -2.46,  -8.78>, < -2.25,   3.81,  -7.67>, Rcell pigment {Black}}
cylinder {<  3.89,  -2.74,  -7.15>, <  3.89,   3.52,  -6.05>, Rcell pigment {Black}}
cylinder {<  2.25,  -3.81,  -1.10>, <  2.25,   2.46,   0.00>, Rcell pigment {Black}}
cylinder {< -3.89,  -3.52,  -2.73>, <  2.25,  -3.81,  -1.10>, Rcell pigment {Black}}
cylinder {< -2.25,  -2.46,  -8.78>, <  3.89,  -2.74,  -7.15>, Rcell pigment {Black}}
cylinder {< -2.25,   3.81,  -7.67>, <  3.89,   3.52,  -6.05>, Rcell pigment {Black}}
cylinder {< -3.89,   2.74,  -1.62>, <  2.25,   2.46,   0.00>, Rcell pigment {Black}}
atom(< -1.95,  -1.76,  -3.56>, 0.56, rgb <0.74, 0.62, 0.89>, 0.0, ase3) // #0
atom(< -1.12,   1.90,  -6.03>, 0.56, rgb <0.74, 0.62, 0.89>, 0.0, ase3) // #1
atom(<  1.95,  -1.37,  -5.77>, 0.56, rgb <0.74, 0.62, 0.89>, 0.0, ase3) // #2
atom(<  1.12,   1.23,  -2.19>, 0.56, rgb <0.74, 0.62, 0.89>, 0.0, ase3) // #3
atom(<  1.12,  -1.90,  -2.75>, 0.56, rgb <0.74, 0.62, 0.89>, 0.0, ase3) // #4
atom(<  1.95,   1.76,  -5.22>, 0.56, rgb <0.74, 0.62, 0.89>, 0.0, ase3) // #5
atom(< -1.95,   1.37,  -3.00>, 0.56, rgb <0.74, 0.62, 0.89>, 0.0, ase3) // #6
atom(< -1.12,  -1.23,  -6.58>, 0.56, rgb <0.74, 0.62, 0.89>, 0.0, ase3) // #7
atom(< -3.89,  -3.52,  -2.73>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #8
atom(< -0.82,  -0.53,  -1.36>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #9
atom(<  0.00,  -3.13,  -4.94>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #10
atom(< -3.07,   0.14,  -5.20>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #11
atom(<  2.25,  -3.81,  -1.10>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #12
atom(<  3.07,  -0.14,  -3.58>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #13
atom(< -3.89,   2.74,  -1.62>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #14
atom(<  0.00,   3.13,  -3.84>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #15
atom(<  2.25,   2.46,   0.00>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #16
atom(< -2.25,  -2.46,  -8.78>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #17
atom(<  0.82,   0.53,  -7.41>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #18
atom(<  3.89,  -2.74,  -7.15>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #19
atom(< -2.25,   3.81,  -7.67>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #20
atom(<  3.89,   3.52,  -6.05>, 0.44, rgb <0.45, 0.68, 0.34>, 0.0, ase3) // #21
cylinder {< -3.89,  -3.52,  -2.73>, < -2.92,  -2.64,  -3.14>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.95,  -1.76,  -3.56>, < -2.92,  -2.64,  -3.14>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -0.82,  -0.53,  -1.36>, < -1.39,  -1.15,  -2.46>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.95,  -1.76,  -3.56>, < -1.39,  -1.15,  -2.46>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -0.82,  -0.53,  -1.36>, <  0.15,   0.35,  -1.78>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.12,   1.23,  -2.19>, <  0.15,   0.35,  -1.78>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -0.82,  -0.53,  -1.36>, <  0.15,  -1.22,  -2.05>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.12,  -1.90,  -2.75>, <  0.15,  -1.22,  -2.05>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -0.82,  -0.53,  -1.36>, < -1.39,   0.42,  -2.18>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.95,   1.37,  -3.00>, < -1.39,   0.42,  -2.18>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -3.13,  -4.94>, < -0.97,  -2.45,  -4.25>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.95,  -1.76,  -3.56>, < -0.97,  -2.45,  -4.25>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -3.13,  -4.94>, <  0.97,  -2.25,  -5.36>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.95,  -1.37,  -5.77>, <  0.97,  -2.25,  -5.36>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -3.13,  -4.94>, <  0.56,  -2.52,  -3.84>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.12,  -1.90,  -2.75>, <  0.56,  -2.52,  -3.84>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,  -3.13,  -4.94>, < -0.56,  -2.18,  -5.76>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.12,  -1.23,  -6.58>, < -0.56,  -2.18,  -5.76>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -3.07,   0.14,  -5.20>, < -2.51,  -0.81,  -4.38>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.95,  -1.76,  -3.56>, < -2.51,  -0.81,  -4.38>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -3.07,   0.14,  -5.20>, < -2.10,   1.02,  -5.61>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.12,   1.90,  -6.03>, < -2.10,   1.02,  -5.61>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -3.07,   0.14,  -5.20>, < -2.51,   0.76,  -4.10>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.95,   1.37,  -3.00>, < -2.51,   0.76,  -4.10>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -3.07,   0.14,  -5.20>, < -2.10,  -0.54,  -5.89>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.12,  -1.23,  -6.58>, < -2.10,  -0.54,  -5.89>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  2.25,  -3.81,  -1.10>, <  1.69,  -2.86,  -1.93>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.12,  -1.90,  -2.75>, <  1.69,  -2.86,  -1.93>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  3.07,  -0.14,  -3.58>, <  2.51,  -0.76,  -4.67>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.95,  -1.37,  -5.77>, <  2.51,  -0.76,  -4.67>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  3.07,  -0.14,  -3.58>, <  2.10,   0.54,  -2.89>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.12,   1.23,  -2.19>, <  2.10,   0.54,  -2.89>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  3.07,  -0.14,  -3.58>, <  2.10,  -1.02,  -3.16>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.12,  -1.90,  -2.75>, <  2.10,  -1.02,  -3.16>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  3.07,  -0.14,  -3.58>, <  2.51,   0.81,  -4.40>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.95,   1.76,  -5.22>, <  2.51,   0.81,  -4.40>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -3.89,   2.74,  -1.62>, < -2.92,   2.06,  -2.31>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.95,   1.37,  -3.00>, < -2.92,   2.06,  -2.31>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   3.13,  -3.84>, < -0.56,   2.52,  -4.93>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.12,   1.90,  -6.03>, < -0.56,   2.52,  -4.93>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   3.13,  -3.84>, <  0.56,   2.18,  -3.01>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.12,   1.23,  -2.19>, <  0.56,   2.18,  -3.01>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   3.13,  -3.84>, <  0.97,   2.45,  -4.53>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.95,   1.76,  -5.22>, <  0.97,   2.45,  -4.53>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.00,   3.13,  -3.84>, < -0.97,   2.25,  -3.42>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.95,   1.37,  -3.00>, < -0.97,   2.25,  -3.42>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  2.25,   2.46,   0.00>, <  1.69,   1.84,  -1.10>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.12,   1.23,  -2.19>, <  1.69,   1.84,  -1.10>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -2.25,  -2.46,  -8.78>, < -1.69,  -1.84,  -7.68>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.12,  -1.23,  -6.58>, < -1.69,  -1.84,  -7.68>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.82,   0.53,  -7.41>, < -0.15,   1.22,  -6.72>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.12,   1.90,  -6.03>, < -0.15,   1.22,  -6.72>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.82,   0.53,  -7.41>, <  1.39,  -0.42,  -6.59>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.95,  -1.37,  -5.77>, <  1.39,  -0.42,  -6.59>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.82,   0.53,  -7.41>, <  1.39,   1.15,  -6.32>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.95,   1.76,  -5.22>, <  1.39,   1.15,  -6.32>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  0.82,   0.53,  -7.41>, < -0.15,  -0.35,  -7.00>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.12,  -1.23,  -6.58>, < -0.15,  -0.35,  -7.00>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  3.89,  -2.74,  -7.15>, <  2.92,  -2.06,  -6.46>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.95,  -1.37,  -5.77>, <  2.92,  -2.06,  -6.46>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {< -2.25,   3.81,  -7.67>, < -1.69,   2.86,  -6.85>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {< -1.12,   1.90,  -6.03>, < -1.69,   2.86,  -6.85>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
cylinder {<  3.89,   3.52,  -6.05>, <  2.92,   2.64,  -5.63>, Rbond texture{pigment {color rgb <0.45, 0.68, 0.34> transmit 0.0} finish{ase3}}}
cylinder {<  1.95,   1.76,  -5.22>, <  2.92,   2.64,  -5.63>, Rbond texture{pigment {color rgb <0.74, 0.62, 0.89> transmit 0.0} finish{ase3}}}
// no constraints
