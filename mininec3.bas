100 rem ****** mininec(3) **********  NOSC code 822 (jcl change 9) 11-26-86
101 rem Original by J.C.Logan & J.W.Rockway @ Navel Ocean Systems Center
102 rem Converted for Chipmunk Basic 3.6.8
103 rem  2019-10-21  by rhn@nicholson.com 
105 print "mininec3 v.05 for Chipmunk Basic"
104 rem
110 defint i,j,k,n
120 dim k(6,2),q(14)
125 mw = 500 : rem increase number of wires and segments
126 mp = 500
130 rem ----- maximum number of segments (pulses + 2 * wires) = 150
140 rem ms = 150
141 ms = 2 * mw + mp
150 dim x(ms),y(ms),z(ms)
160 rem ----- maximum number of wires = 50
170 rem mw = 50
180 dim a(mw),ca(mw),cb(mw),cg(mw),j1(mw),j2(mw,2),n(mw,2),s(mw)
190 rem ----- maximum number of loads = 11
200 ml = 11
210 rem ----- maximum order of s-parameter loads = 8
220 ma = 8
230 dim la(2,11,8),lp(11),ls(11)
240 rem ----- maximum number of media = 6
250 mm = 6
260 rem ----- h must be dimensioned at least 6
270 dim h(6),t(6),u(6),v(6),z1(6),z2(6)
280 rem ----- maximum number of pulses = 50 
290 rem mp = 50
300 dim c%(mp,2),ci(mp),cr(mp),p(mp),w%(mp)
310 dim zr(mp,mp),zi(mp,mp)
320 rem ---- arrays e,l & m dimensioned to mw+mp=100
330 dim e(mw+mp),l(mw+mp),m(mw+mp)
340 rem graphics color 2,0
350 goto 15060
360 rem ********** kernel evaluation of integrals i2 & i3 **********
370 if k < 0 then 420
380 x3 = x2+t*(v1-x2)
390 y3 = y2+t*(v2-y2)
400 z3 = z2+t*(v3-z2)
410 goto 450
420 x3 = v1+t*(x2-v1)
430 y3 = v2+t*(y2-v2)
440 z3 = v3+t*(z2-v3)
450 d3 = x3*x3+y3*y3+z3*z3
460 rem ----- mod for small radius to wavelength ratio
470 if a(p4) <= srm then d = sqr(d3) : goto 580
480 d = d3+a2
490 if d > 0 then d = sqr(d)
500 rem ----- criteria for using reduced kernel
510 if i6 = 0 then 580
520 rem ----- exact kernel calculation with elliptic integral
530 b = d3/(d3+4*a2)
540 w0 = c0+b*(c1+b*(c2+b*(c3+b*c4)))
550 w1 = c5+b*(c6+b*(c7+b*(c8+b*c9)))
560 v0 = (w0-w1*log(b))*sqr(1-b)
570 t3 = t3+(v0+log(d3/(64*a2))/2)/p/a(p4)-1/d
580 b1 = d*w
590 rem ----- exp(-j*k*r)/r
600 t3 = t3+cos(b1)/d
610 t4 = t4-sin(b1)/d
620 return
630 rem ***** psi(p1,p2,p3) = t1 + j * t2 **********
640 rem ----- entries required for near field calculation
650 x1 = x0+p1*t5/2
660 y1 = y0+p1*t6/2
670 z1 = z0+p1*t7/2
680 x2 = x1-x(p2)
690 y2 = y1-y(p2)
700 z2 = z1-k*z(p2)
710 v1 = x1-x(p3)
720 v2 = y1-y(p3)
730 v3 = z1-k*z(p3)
740 goto 1440
750 i4 = int(p2)
760 i5 = i4+1
770 x2 = x0-(x(i4)+x(i5))/2
780 y2 = y0-(y(i4)+y(i5))/2
790 z2 = z0-k*(z(i4)+z(i5))/2
800 v1 = x0-x(p3)
810 v2 = y0-y(p3)
820 v3 = z0-k*z(p3)
830 goto 1440
840 x2 = x0-x(p2)
850 y2 = y0-y(p2)
860 z2 = z0-k*z(p2)
870 i4 = int(p3)
880 i5 = i4+1
890 v1 = x0-(x(i4)+x(i5))/2
900 v2 = y0-(y(i4)+y(i5))/2
910 v3 = z0-k*(z(i4)+z(i5))/2
920 goto 1440
930 rem ----- entries required for impedance matrix calculation
940 rem ----- s(m) goes in (x1,y1,z1) for scalar potential
950 rem ----- mod for small radius to wave length ratio
960 fvs = 1
970 if k < 1 then 1030
980 if a(p4) > srm then 1030
990 if (p3 = p2+1 and p1 = (p2+p3)/2) then goto 1000 : else goto 1030
1000 t1 = 2*log(s(p4)/a(p4))
1010 t2 = -w*s(p4)
1020 return
1030 i4 = int(p1)
1040 i5 = i4+1
1050 x1 = (x(i4)+x(i5))/2
1060 y1 = (y(i4)+y(i5))/2
1070 z1 = (z(i4)+z(i5))/2
1080 goto 1220
1090 rem ----- s(m) goes in (x1,y1,z1) for vector potential
1100 rem ----- mod for small radius to wave length ratio
1110 fvs = 0
1120 if k < 1 then 1180
1130 if a(p4) >= srm then 1180
1140 if (i = j and p3 = p2+0.5) then 1150 : else goto 1180
1150 t1 = log(s(p4)/a(p4))
1160 t2 = -w*s(p4)/2
1170 return
1180 x1 = x(p1)
1190 y1 = y(p1)
1200 z1 = z(p1)
1210 rem ----- s(u)-s(m) goes in (x2,y2,z2)
1220 i4 = int(p2)
1230 if i4 = p2 then 1290
1240 i5 = i4+1
1250 x2 = (x(i4)+x(i5))/2-x1
1260 y2 = (y(i4)+y(i5))/2-y1
1270 z2 = k*(z(i4)+z(i5))/2-z1
1280 goto 1330
1290 x2 = x(p2)-x1
1300 y2 = y(p2)-y1
1310 z2 = k*z(p2)-z1
1320 rem ----- s(v)-s(m) goes in (v1,v2,v3)
1330 i4 = int(p3)
1340 if i4 = p3 then 1400
1350 i5 = i4+1
1360 v1 = (x(i4)+x(i5))/2-x1
1370 v2 = (y(i4)+y(i5))/2-y1
1380 v3 = k*(z(i4)+z(i5))/2-z1
1390 goto 1440
1400 v1 = x(p3)-x1
1410 v2 = y(p3)-y1
1420 v3 = k*z(p3)-z1
1430 rem ----- magnitude of s(u) - s(m)
1440 d0 = x2*x2+y2*y2+z2*z2
1450 rem ----- magnitude of s(v) - s(m)
1460 if d0 > 0 then d0 = sqr(d0)
1470 d3 = v1*v1+v2*v2+v3*v3
1480 if d3 > 0 then d3 = sqr(d3)
1490 rem ----- square of wire radius
1500 a2 = a(p4)*a(p4)
1510 rem ----- magnitude of s(v) - s(u)
1520 s4 = (p3-p2)*s(p4)
1530 rem ----- order of integration
1540 rem ----- lth order gaussian quadrature
1550 t1 = 0
1560 t2 = 0
1570 i6 = 0
1580 f2 = 1
1590 l = 7
1600 t = (d0+d3)/s(p4)
1610 rem ----- criteria for exact kernel
1620 if t > 1.1 then 1740
1630 if c$ = "n" then 1740
1640 if j2(w%(i),1) = j2(w%(j),1) then 1690
1650 if j2(w%(i),1) = j2(w%(j),2) then 1690
1660 if j2(w%(i),2) = j2(w%(j),1) then 1690
1670 if j2(w%(i),2) = j2(w%(j),2) then 1690
1680 goto 1740
1690 if a(p4) > srm then 1710
1700 if fvs = 1 then 1000 : else goto 1150
1710 f2 = 2*(p3-p2)
1720 i6 = (1-log(s4/f2/8/a(p4)))/p/a(p4)
1730 goto 1760
1740 if t > 6 then l = 3
1750 if t > 10 then l = 1
1760 i5 = l+l
1770 t3 = 0
1780 t4 = 0
1790 t = (q(l)+0.5)/f2
1800 gosub 370
1810 t = (0.5-q(l))/f2
1820 gosub 370
1830 l = l+1
1840 t1 = t1+q(l)*t3
1850 t2 = t2+q(l)*t4
1860 l = l+1
1870 if l < i5 then 1770
1880 t1 = s4*(t1+i6)
1890 t2 = s4*t2
1900 return
1910 rem ********** complex square root **********
1920 rem ----- w6+i*w7=sqr(z6+i*z7)
1930 t6 = sqr((abs(z6)+sqr(z6*z6+z7*z7))/2)
1940 t7 = abs(z7)/2/t6
1950 if z6 < 0 then 2000
1960 w6 = t6
1970 w7 = t7
1980 if z7 < 0 then w7 = -t7
1990 return
2000 w6 = t7
2010 w7 = t6
2020 if z7 < 0 then w7 = -t6
2030 return
2040 rem ********** impedance matrix calculation **********
2050 if flg = 1 then 4370
2060 if flg = 2 then 4860
2070 rem ----- begin matrix fill time calculation
2080 ot$ = time$
2090 q$ = "matrix fill  "
2100 print
2110 print "begin ";q$
2120 rem ----- zero impedance matrix
2130 for i = 1 to n
2140 for j = 1 to n
2150 zr(i,j) = 0
2160 zi(i,j) = 0
2170 next j
2180 next i
2190 rem ----- compute row i of matrix (observation loop)
2200 for i = 1 to n
2210 i1 = abs(c%(i,1))
2220 i2 = abs(c%(i,2))
2230 f4 = sgn(c%(i,1))*s(i1)
2240 f5 = sgn(c%(i,2))*s(i2)
2250 rem ----- r(m + 1/2) - r(m - 1/2) has components (t5,t6,t7)
2260 t5 = f4*ca(i1)+f5*ca(i2)
2270 t6 = f4*cb(i1)+f5*cb(i2)
2280 t7 = f4*cg(i1)+f5*cg(i2)
2290 if c%(i,1) = -c%(i,2) then t7 = s(i1)*(cg(i1)+cg(i2))
2300 rem ----- compute column j of row i (source loop)
2310 for j = 1 to n
2320 j1 = abs(c%(j,1))
2330 j2 = abs(c%(j,2))
2340 f4 = sgn(c%(j,1))
2350 f5 = sgn(c%(j,2))
2360 f6 = 1
2370 f7 = 1
2380 rem ----- image loop
2390 for k = 1 to g step -2
2400 if c%(j,1) <> -c%(j,2) then 2440
2410 if k < 0 then 3410
2420 f6 = f4
2430 f7 = f5
2440 f8 = 0
2450 if k < 0 then 2570
2460 rem ----- set flag to avoid redunant calculations
2470 if i1 <> i2 then 2550
2480 if (ca(i1)+cb(i1)) = 0 then 2500
2490 if c%(i,1) <> c%(i,2) then 2550
2500 if j1 <> j2 then 2550
2510 if (ca(j1)+cb(j1)) = 0 then 2530
2520 if c%(j,1) <> c%(j,2) then 2550
2530 if i1 = j1 then f8 = 1
2540 if i = j then f8 = 2
2550 if zr(i,j) <> 0 then 3260
2560 rem ----- compute psi(m,n,n+1/2)
2570 p1 = 2*w%(i)+i-1
2580 p2 = 2*w%(j)+j-1
2590 p3 = p2+0.5
2600 p4 = j2
2610 gosub 1110
2620 u1 = f5*t1
2630 u2 = f5*t2
2640 rem ----- compute psi(m,n-1/2,n)
2650 p3 = p2
2660 p2 = p2-0.5
2670 p4 = j1
2680 if f8 < 2 then gosub 1110
2690 v1 = f4*t1
2700 v2 = f4*t2
2710 rem ----- s(n+1/2)*psi(m,n,n+1/2) + s(n-1/2)*psi(m,n-1/2,n)
2720 x3 = u1*ca(j2)+v1*ca(j1)
2730 y3 = u1*cb(j2)+v1*cb(j1)
2740 z3 = (f7*u1*cg(j2)+f6*v1*cg(j1))*k
2750 rem ----- real part of vector potential contribution
2760 d1 = w2*(x3*t5+y3*t6+z3*t7)
2770 x3 = u2*ca(j2)+v2*ca(j1)
2780 y3 = u2*cb(j2)+v2*cb(j1)
2790 z3 = (f7*u2*cg(j2)+f6*v2*cg(j1))*k
2800 rem ----- imaginary part of vector potential contribution
2810 d2 = w2*(x3*t5+y3*t6+z3*t7)
2820 rem ----- compute psi(m+1/2,n,n+1)
2830 p1 = p1+0.5
2840 if f8 = 2 then p1 = p1-1
2850 p2 = p3
2860 p3 = p3+1
2870 p4 = j2
2880 if f8 <> 1 then 2920
2890 u5 = f5*u1+t1
2900 u6 = f5*u2+t2
2910 goto 3000
2920 gosub 960
2930 if f8 < 2 then 2970
2940 u1 = (2*t1-4*u1*f5)/s(j1)
2950 u2 = (2*t2-4*u2*f5)/s(j1)
2960 goto 3230
2970 u5 = t1
2980 u6 = t2
2990 rem ----- compute psi(m-1/2,n,n+1)
3000 p1 = p1-1
3010 gosub 960
3020 u1 = (t1-u5)/s(j2)
3030 u2 = (t2-u6)/s(j2)
3040 rem ----- compute psi(m+1/2,n-1,n)
3050 p1 = p1+1
3060 p3 = p2
3070 p2 = p2-1
3080 p4 = j1
3090 gosub 960
3100 u3 = t1
3110 u4 = t2
3120 rem ----- compute psi(m-1/2,n-1,n)
3130 if f8 < 1 then 3170
3140 t1 = u5
3150 t2 = u6
3160 goto 3200
3170 p1 = p1-1
3180 gosub 960
3190 rem ----- gradient of scalar potential contribution
3200 u1 = u1+(u3-t1)/s(j1)
3210 u2 = u2+(u4-t2)/s(j1)
3220 rem ----- sum into impedance matrix
3230 zr(i,j) = zr(i,j)+k*(d1+u1)
3240 zi(i,j) = zi(i,j)+k*(d2+u2)
3250 rem ----- avoid redunant calculations
3260 if j < i then 3410
3270 if f8 = 0 then 3410
3280 zr(j,i) = zr(i,j)
3290 zi(j,i) = zi(i,j)
3300 rem ----- segments on same wire same distance apart have same z
3310 p1 = j+1
3320 if p1 > n then 3410
3330 if c%(p1,1) <> c%(p1,2) then 3410
3340 if c%(p1,2) = c%(j,2) then 3370
3350 if c%(p1,2) <> -c%(j,2) then 3410
3360 if (ca(j2)+cb(j2)) <> 0 then 3410
3370 p2 = i+1
3380 if p2 > n then 3410
3390 zr(p2,p1) = zr(i,j)
3400 zi(p2,p1) = zi(i,j)
3410 next k
3420 next j
3430 pct = i/n
3440 gosub 16080
3450 next i
3460 rem ----- end matrix fill time calculation
3470 t$ = time$
3480 gosub 15980
3490 print #3," "
3500 print #3,"fill matrix  : ";t$
3510 rem ********** addition of loads **********
3520 if nl = 0 then 3860
3530 f5 = 2*p*f
3540 for i = 1 to nl
3550 if l$ = "n" then 3750
3560 rem ----- s-parameter loads
3570 u1 = 0
3580 u2 = 0
3590 d1 = 0
3600 d2 = 0
3610 s = 1
3620 for j = 0 to ls(i) step 2
3630 u1 = u1+la(1,i,j)*s*f5^j
3640 d1 = d1+la(2,i,j)*s*f5^j
3650 l = j+1
3660 u2 = u2+la(1,i,l)*s*f5^l
3670 d2 = d2+la(2,i,l)*s*f5^l
3680 s = -s
3690 next j
3700 j = lp(i)
3710 d = d1*d1+d2*d2
3720 li = (u2*d1-d2*u1)/d
3730 lr = (u1*d1+u2*d2)/d
3740 goto 3780
3750 lr = la(1,i,1)
3760 li = la(2,i,1)
3770 j = lp(i)
3780 f2 = 1/m
3790 if c%(j,1) <> -c%(j,2) then 3810
3800 if k < 0 then f2 = 2/m
3810 zr(j,j) = zr(j,j)+f2*li
3820 zi(j,j) = zi(j,j)-f2*lr
3830 next i
3840 rem ********** impedance matrix factorization **********
3850 rem ----- begin matrix factor time calculation
3860 ot$ = time$
3870 q$ = "factor matrix"
3880 print
3890 print "begin ";q$;
3900 x = n
3910 pctn = x*(x-1)*(x+x-1)
3920 for k = 1 to n-1
3930 rem ----- search for pivot
3940 t = zr(k,k)*zr(k,k)+zi(k,k)*zi(k,k)
3950 i1 = k
3960 for i = k+1 to n
3970 t1 = zr(i,k)*zr(i,k)+zi(i,k)*zi(i,k)
3980 if t1 < t then 4010
3990 i1 = i
4000 t = t1
4010 next i
4020 rem ----- exchange rows k and i1
4030 if i1 = k then 4120
4040 for j = 1 to n
4050 t1 = zr(k,j)
4060 t2 = zi(k,j)
4070 zr(k,j) = zr(i1,j)
4080 zi(k,j) = zi(i1,j)
4090 zr(i1,j) = t1
4100 zi(i1,j) = t2
4110 next j
4120 p(k) = i1
4130 rem ----- subtract row k from rows k+1 to n
4140 for i = k+1 to n
4150 rem ----- compute multiplier l(i,k)
4160 t1 = (zr(i,k)*zr(k,k)+zi(i,k)*zi(k,k))/t
4170 t2 = (zi(i,k)*zr(k,k)-zr(i,k)*zi(k,k))/t
4180 zr(i,k) = t1
4190 zi(i,k) = t2
4200 rem ----- subtract row k from row i
4210 for j = k+1 to n
4220 zr(i,j) = zr(i,j)-(zr(k,j)*t1-zi(k,j)*t2)
4230 zi(i,j) = zi(i,j)-(zr(k,j)*t2+zi(k,j)*t1)
4240 next j
4250 next i
4260 x = n-k
4270 pct = 1-x*(x-1)*(x+x-1)/pctn
4280 gosub 16080
4290 next k
4300 rem ----- end matrix factor time calculation
4310 t$ = time$
4320 gosub 15980
4330 print
4340 print #3,"factor matrix: ";t$
4350 rem ********** solve **********
4360 rem ----- compute right hand side
4370 for i = 1 to n
4380 cr(i) = 0
4390 ci(i) = 0
4400 next i
4410 for j = 1 to ns
4420 f2 = 1/m
4430 if c%(e(j),1) = -c%(e(j),2) then f2 = 2/m
4440 cr(e(j)) = f2*m(j)
4450 ci(e(j)) = -f2*l(j)
4460 next j
4470 rem ----- permute excitation
4480 for k = 1 to n-1
4490 i1 = p(k)
4500 if i1 = k then 4570
4510 t1 = cr(k)
4520 t2 = ci(k)
4530 cr(k) = cr(i1)
4540 ci(k) = ci(i1)
4550 cr(i1) = t1
4560 ci(i1) = t2
4570 next k
4580 rem ----- forward elimination
4590 for i = 2 to n
4600 t1 = 0
4610 t2 = 0
4620 for j = 1 to i-1
4630 t1 = t1+zr(i,j)*cr(j)-zi(i,j)*ci(j)
4640 t2 = t2+zr(i,j)*ci(j)+zi(i,j)*cr(j)
4650 next j
4660 cr(i) = cr(i)-t1
4670 ci(i) = ci(i)-t2
4680 next i
4690 rem ----- back substitution
4700 for i = n to 1 step -1
4710 t1 = 0
4720 t2 = 0
4730 if i = n then 4780
4740 for j = i+1 to n
4750 t1 = t1+zr(i,j)*cr(j)-zi(i,j)*ci(j)
4760 t2 = t2+zr(i,j)*ci(j)+zi(i,j)*cr(j)
4770 next j
4780 t = zr(i,i)*zr(i,i)+zi(i,i)*zi(i,i)
4790 t1 = cr(i)-t1
4800 t2 = ci(i)-t2
4810 cr(i) = (t1*zr(i,i)+t2*zi(i,i))/t
4820 ci(i) = (t2*zr(i,i)-t1*zi(i,i))/t
4830 next i
4840 flg = 2
4850 rem ********** source data **********
4860 print #3," "
4870 print #3,b$;"    source data     ";b$
4880 pwr = 0
4890 for i = 1 to ns
4900 cr = cr(e(i))
4910 ci = ci(e(i))
4920 t = cr*cr+ci*ci
4930 t1 = (l(i)*cr+m(i)*ci)/t
4940 t2 = (m(i)*cr-l(i)*ci)/t
4950 o2 = (l(i)*cr+m(i)*ci)/2
4960 pwr = pwr+o2
4970 print #3,"pulse ";e(i),"voltage = (";l(i);",";m(i);"j)"
4980 print #3," ","current = (";cr;",";ci;"j)"
4990 print #3," ","impedance = (";t1;",";t2;"j)"
5000 print #3," ","power = ";o2;" watts"
5010 next i
5020 if ns > 1 then print #3," "
5030 if ns > 1 then print #3,"total power = ";pwr;"watts"
5040 return
5050 rem ********** print currents **********
5060 gosub 2050
5070 s$ = "n"
5080 print #3," "
5090 print #3,b$;"    current data    ";b$
5100 for k = 1 to nw
5110 if s$ = "y" then 5160
5120 print #3," "
5130 print #3,"wire no. ";k;":"
5140 print #3,"pulse","real","imaginary","magnitude","phase"
5150 print #3," no.","(amps)","(amps)","(amps)","(degrees)"
5160 n1 = n(k,1)
5170 n2 = n(k,2)
5180 i = n1
5190 c = c%(i,1)
5200 if (n1 = 0 and n2 = 0) then c = k
5210 if g = 1 then 5240
5220 if (j1(k) = -1 and n1 > n2) then n2 = n1
5230 if j1(k) = -1 then 5340
5240 e% = 1
5250 gosub 5810
5260 i2 = i1
5270 j2 = j1
5280 gosub 6160
5290 if s$ = "n" then print #3,i$,i1;tab (29);j1;tab (43);s1;tab (57);s2
5300 if s$ = "y" then print #1,i1;",";j1;",";s1;",";s2
5310 if n1 = 0 then 5410
5320 if c = k then 5340
5330 if i$ = "j" then n1 = n1+1
5340 for i = n1 to n2-1
5350 i2 = cr(i)
5360 j2 = ci(i)
5370 gosub 6160
5380 if s$ = "n" then print #3,i,cr(i);tab (29);ci(i);tab (43);s1;tab (57);s2
5390 if s$ = "y" then print #1,cr(i);",";ci(i);",";s1;",";s2
5400 next i
5410 i = n2
5420 c = c%(i,2)
5430 if (n1 = 0 and n2 = 0) then c = k
5440 if g = 1 then 5460
5450 if j1(k) = 1 then 5520
5460 e% = 2
5470 gosub 5810
5480 if (n1 = 0 and n2 = 0) then 5580
5490 if n1 > n2 then 5580
5500 if c = k then 5520
5510 if i$ = "j" then 5580
5520 i2 = cr(n2)
5530 j2 = ci(n2)
5540 gosub 6160
5550 if s$ = "n" then print #3,n2,cr(n2);tab (29);ci(n2);tab (43);s1;tab (57);s2
5560 if s$ = "y" then print #1,cr(n2);",";ci(n2);",";s1;",";s2
5570 if j1(k) = 1 then 5630
5580 i2 = i1
5590 j2 = j1
5600 gosub 6160
5610 if s$ = "n" then print #3,i$,i1;tab (29);j1;tab (43);s1;tab (57);s2
5620 if s$ = "y" then print #1,i1;",";j1;",";s1;",";s2
5630 if s$ = "y" then print #1," 1 , 1 , 1 , 1"
5640 next k
5650 if s$ = "y" then 5780
5660 print
5670 input "save currents to a file (y/n) > ";s$
5680 if s$ = "n" then 5790
5690 if s$ <> "y" then 5660
5700 print #3," "
5710 input "filename (name.out) > ";f$
5720 if left$(right$(f$,4),1) = "." then 5730 : else f$ = f$+".out"
5730 if o$ > "c" then print #3,"filename (name.out): ";f$
5740 open f$ for output as #1
5750 print #3," "
5760 print #1,nw;",";pwr;",c"
5770 goto 5100
5780 close #1
5790 return
5800 rem ----- sort junction currents
5810 i$ = "e"
5820 i1 = 0
5830 j1 = 0
5840 if (c = k or c = 0) then goto 5890
5850 i$ = "j"
5860 i1 = cr(i)
5870 j1 = ci(i)
5880 rem ----- check for other overlapping wires
5890 for j = 1 to nw
5900 if j = k then goto 6130
5910 l1 = n(j,1)
5920 l2 = n(j,2)
5930 if e% = 2 then goto 5990
5940 co = c%(l1,1)
5950 ct = c%(l2,2)
5960 l3 = l1
5970 l4 = l2
5980 goto 6030
5990 co = c%(l2,2)
6000 ct = c%(l1,1)
6010 l3 = l2
6020 l4 = l1
6030 if co = -k then 6050
6040 goto 6080
6050 i1 = i1-cr(l3)
6060 j1 = j1-ci(l3)
6070 i$ = "j"
6080 if ct = k then 6100
6090 goto 6130
6100 i1 = i1+cr(l4)
6110 j1 = j1+ci(l4)
6120 i$ = "j"
6130 next j
6140 return
6150 rem ----- calculate s1 and s2
6160 i3 = i2*i2
6170 j3 = j2*j2
6180 if (i3 > 0 or j3 > 0) then 6210
6190 s1 = 0
6200 goto 6220
6210 s1 = sqr(i3+j3)
6220 if i2 <> 0 then 6250
6230 s2 = 0
6240 return
6250 s2 = arctan(j2/i2)/p0
6260 if i2 > 0 then return
6270 s2 = s2+sgn(j2)*180
6280 return
6290 rem ********** far field calculation **********
6300 if flg < 2 then gosub 2050
6310 o2 = pwr
6320 rem ----- tabulate impedance
6330 if nm = 0 then 6430
6340 for i = 1 to nm
6350 z6 = t(i)
6360 z7 = -v(i)/(2*p*f*8.850000E-06)
6370 rem ----- form impedance=1/sqr(dielectric constant)
6380 gosub 1930
6390 d = w6*w6+w7*w7
6400 z1(i) = w6/d
6410 z2(i) = -w7/d
6420 next i
6430 print #3," "
6440 print #3,b$;"     far field      ";b$
6450 print #3," "
6460 rem ----- input variables for far field calculation
6470 input "calculate pattern in dbi or volts/meter (d/v) > ";p$
6480 if p$ = "d" then 6640
6490 if p$ <> "v" then 6470
6500 f1 = 1
6510 print
6520 print "present power level =  ";pwr;" watts"
6530 input "change power level (y/n) > ";a$
6540 if a$ = "n" then 6590
6550 if a$ <> "y" then 6530
6560 input "new power level (watts) > ";o2
6570 if o$ > "c" then print #3,"new power level = ";o2
6580 goto 6530
6590 if (o2 < 0 or o2 = 0) then o2 = pwr
6600 f1 = sqr(o2/pwr)
6610 print
6620 input "radial distance (meters) > ";rd
6630 if rd < 0 then rd = 0
6640 a$ = "zenith angle : initial,increment,number > "
6650 print a$;
6660 input za,zc,nz
6670 if nz = 0 then nz = 1
6680 if o$ > "c" then print #3,a$;": ";za;",";zc;",";nz
6690 a$ = "azimuth angle: initial,increment,number > "
6700 print a$;
6710 input aa,ac,na
6720 if na = 0 then na = 1
6730 if o$ > "c" then print #3,a$;": ";aa;",";ac;",";na
6740 print #3," "
6750 rem ********** file far field data **********
6760 input "file pattern (y/n) > ";s$
6770 if s$ = "n" then 6850
6780 if s$ <> "y" then 6760
6790 print #3," "
6800 input "filename (name.out) > ";f$
6810 if left$(right$(f$,4),1) = "." then 6820 : else f$ = f$+".out"
6820 if o$ > "c" then print #3,"filename (name.out): ";f$
6830 open f$ for output as #1
6840 print #1,na*nz;",";o2;",";p$
6850 print #3," "
6860 k9 = 0.016678/pwr
6870 rem ----- pattern header
6880 print #3,b$;"    pattern data    ";b$
6890 if p$ = "v" then goto 6940
6900 print #3,"zenith","azimuth","vertical","horizontal","total"
6910 a$ = "pattern (db)"
6920 print #3," angle"," angle",a$,a$,a$
6930 goto 7010
6940 if rd > 0 then print #3,tab (15);"radial distance = ";rd;" meters"
6950 print #3,tab (15);"power level = ";pwr*f1*f1;" watts"
6960 print #3,"zenith   azimuth","     e(theta)     ","     e(phi)"
6970 a$ = " mag(v/m)    phase(deg)"
6980 print #3," angle    angle",a$,a$
6990 if s$ = "y" then print #1,rd
7000 rem ----- loop over azimuth angle
7010 q1 = aa
7020 for i1 = 1 to na
7030 u3 = q1*p0
7040 v1 = -sin(u3)
7050 v2 = cos(u3)
7060 rem ----- loop over zenith angle
7070 q2 = za
7080 for i2 = 1 to nz
7090 u4 = q2*p0
7100 r3 = cos(u4)
7110 t3 = -sin(u4)
7120 t1 = r3*v2
7130 t2 = -r3*v1
7140 r1 = -t3*v2
7150 r2 = t3*v1
7160 x1 = 0
7170 y1 = 0
7180 z1 = 0
7190 x2 = 0
7200 y2 = 0
7210 z2 = 0
7220 rem ----- image loop
7230 for k = 1 to g step -2
7240 for i = 1 to n
7250 if k > 0 then 7270
7260 if c%(i,1) = -c%(i,2) then 8220
7270 j = 2*w%(i)-1+i
7280 rem ----- for each end of pulse compute a contribution to e-field
7290 for f5 = 1 to 2
7300 l = abs(c%(i,f5))
7310 f3 = sgn(c%(i,f5))*w*s(l)/2
7320 if c%(i,1) <> -c%(i,2) then 7340
7330 if f3 < 0 then 8210
7340 if k = 1 then 7370
7350 if nm <> 0 then 7560
7360 rem ----- standard case
7370 s2 = w*(x(j)*r1+y(j)*r2+z(j)*k*r3)
7380 s1 = cos(s2)
7390 s2 = sin(s2)
7400 b1 = f3*(s1*cr(i)-s2*ci(i))
7410 b2 = f3*(s1*ci(i)+s2*cr(i))
7420 if c%(i,1) = -c%(i,2) then 7510
7430 x1 = x1+k*b1*ca(l)
7440 x2 = x2+k*b2*ca(l)
7450 y1 = y1+k*b1*cb(l)
7460 y2 = y2+k*b2*cb(l)
7470 z1 = z1+b1*cg(l)
7480 z2 = z2+b2*cg(l)
7490 goto 8210
7500 rem ----- grounded ends
7510 z1 = z1+2*b1*cg(l)
7520 z2 = z2+2*b2*cg(l)
7530 goto 8210
7540 rem ----- real ground case
7550 rem ----- begin by finding specular distance
7560 t4 = 100000
7570 if r3 = 0 then 7590
7580 t4 = -z(j)*t3/r3
7590 b9 = t4*v2+x(j)
7600 if tb = 1 then 7640
7610 b9 = b9*b9+(y(j)-t4*v1)^2
7620 if b9 > 0 then b9 = sqr(b9) : else goto 7640
7630 rem ----- search for the corresponding medium
7640 j2 = nm
7650 for j1 = nm to 1 step -1
7660 if b9 > u(j1) then 7680
7670 j2 = j1
7680 next j1
7690 rem ----- obtain impedance at specular point
7700 z4 = z1(j2)
7710 z5 = z2(j2)
7720 rem ----- if present include ground screen impedance in parallel
7730 if nr = 0 then 7850
7740 if b9 > u(1) then 7850
7750 r = b9+nr*rr
7760 z8 = w*r*log(r/(nr*rr))/nr
7770 s8 = -z5*z8
7780 s9 = z4*z8
7790 t8 = z4
7800 t9 = z5+z8
7810 d = t8*t8+t9*t9
7820 z4 = (s8*t8+s9*t9)/d
7830 z5 = (s9*t8-s8*t9)/d
7840 rem ----- form sqr(1-z^2*sin^2)
7850 z6 = 1-(z4*z4-z5*z5)*t3*t3
7860 z7 = -(2*z4*z5)*t3*t3
7870 gosub 1930
7880 rem ----- vertical reflection coefficient
7890 s8 = r3-(w6*z4-w7*z5)
7900 s9 = -(w6*z5+w7*z4)
7910 t8 = r3+(w6*z4-w7*z5)
7920 t9 = w6*z5+w7*z4
7930 d = t8*t8+t9*t9
7940 v8 = (s8*t8+s9*t9)/d
7950 v9 = (s9*t8-s8*t9)/d
7960 rem ----- horizontal reflection coefficient
7970 s8 = w6-r3*z4
7980 s9 = w7-r3*z5
7990 t8 = w6+r3*z4
8000 t9 = w7+r3*z5
8010 d = t8*t8+t9*t9
8020 h8 = (s8*t8+s9*t9)/d-v8
8030 h9 = (s9*t8-s8*t9)/d-v9
8040 rem ----- compute contribution to sum
8050 s2 = w*(x(j)*r1+y(j)*r2-(z(j)-2*h(j2))*r3)
8060 s1 = cos(s2)
8070 s2 = sin(s2)
8080 b1 = f3*(s1*cr(i)-s2*ci(i))
8090 b2 = f3*(s1*ci(i)+s2*cr(i))
8100 w6 = b1*v8-b2*v9
8110 w7 = b1*v9+b2*v8
8120 d = ca(l)*v1+cb(l)*v2
8130 z6 = d*(b1*h8-b2*h9)
8140 z7 = d*(b1*h9+b2*h8)
8150 x1 = x1-(ca(l)*w6+v1*z6)
8160 x2 = x2-(ca(l)*w7+v1*z7)
8170 y1 = y1-(cb(l)*w6+v2*z6)
8180 y2 = y2-(cb(l)*w7+v2*z7)
8190 z1 = z1+cg(l)*w6
8200 z2 = z2+cg(l)*w7
8210 next f5
8220 next i
8230 next k
8240 h2 = -(x1*t1+y1*t2+z1*t3)*g0
8250 h1 = (x2*t1+y2*t2+z2*t3)*g0
8260 x4 = -(x1*v1+y1*v2)*g0
8270 x3 = (x2*v1+y2*v2)*g0
8280 if p$ = "d" then 8360
8290 if rd = 0 then 8510
8300 h1 = h1/rd
8310 h2 = h2/rd
8320 x3 = x3/rd
8330 x4 = x4/rd
8340 goto 8510
8350 rem ----- pattern in db
8360 p1 = -999
8370 p2 = p1
8380 p3 = p1
8390 t1 = k9*(h1*h1+h2*h2)
8400 t2 = k9*(x3*x3+x4*x4)
8410 t3 = t1+t2
8420 rem ----- calculate values in db
8430 if t1 > 1.000000E-30 then p1 = 4.343*log(t1)
8440 if t2 > 1.000000E-30 then p2 = 4.343*log(t2)
8450 if t3 > 1.000000E-30 then p3 = 4.343*log(t3)
8460 print #3,q2;tab (15);q1;tab (29);p1;tab (43);p2;tab (57);p3
8470 if s$ = "y" then print #1,q2;",";q1;",";p1;",";p2;",";p3
8480 goto 8750
8490 rem ----- pattern in volts/meter
8500 rem ----- magnitude and phase of e(theta)
8510 s1 = 0
8520 if (h1 = 0 and h2 = 0) then 8540
8530 s1 = sqr(h1*h1+h2*h2)
8540 if h1 <> 0 then 8570
8550 s2 = 0
8560 goto 8600
8570 s2 = arctan(h2/h1)/p0
8580 if h1 < 0 then s2 = s2+sgn(h2)*180
8590 rem ----- magnitude and phase of e(phi)
8600 s3 = 0
8610 if (x3 = 0 and x4 = 0) then 8630
8620 s3 = sqr(x3*x3+x4*x4)
8630 if x3 <> 0 then 8660
8640 s4 = 0
8650 goto 8680
8660 s4 = arctan(x4/x3)/p0
8670 if x3 < 0 then s4 = s4+sgn(x4)*180
8680 print #3,using "###.##    ";q2,q1;
8690 print #3,using "       ##.###^^^^";s1*f1;
8700 print #3,using "   ###.##   ";s2;
8710 print #3,using "       ##.###^^^^";s3*f1;
8720 print #3,using "   ###.##";s4
8730 if s$ = "y" then print #1,q2;",";q1;",";s1*f1;",";s2;",";s3*f1;"," s4
8740 rem ----- increment zenith angle
8750 q2 = q2+zc
8760 next i2
8770 rem ----- increment azimuth angle
8780 q1 = q1+ac
8790 next i1
8800 close #1
8810 return
8820 rem ********** near field calculation **********
8830 rem ----- ensure currents have been calculated
8840 if flg < 2 then gosub 2050
8850 o2 = pwr
8860 print #3," "
8870 print #3,b$;"    near fields     ";b$
8880 print #3," "
8890 input "electric or magnetic near fields (e/h) > ";n$
8900 if (n$ = "h" or n$ = "e")goto 8920
8910 goto 8890
8920 print
8930 rem ----- input variables for near field calculation
8940 print "field location(s):"
8950 a$ = "-coordinate (m): initial,increment,number "
8960 print "   x";a$;" > ";
8970 input xi,xc,nx
8980 if nx = 0 then nx = 1
8990 if o$ > "c" then print #3,"x";a$;": ";xi;",";xc;",";nx
9000 print "   y";a$;" > ";
9010 input yi,yc,ny
9020 if ny = 0 then ny = 1
9030 if o$ > "c" then print #3,"y";a$;": ";yi;",";yc;",";ny
9040 print "   z";a$;" > ";
9050 input zi,zc,nz
9060 if nz = 0 then nz = 1
9070 if o$ > "c" then print #3,"z";a$;": ";zi;",";zc;",";nz
9080 f1 = 1
9090 print
9100 print "present power level is ";pwr;" watts"
9110 input "change power level (y/n) > ";a$
9120 if a$ = "n" then 9170
9130 if a$ <> "y" then 9110
9140 input "new power level (watts) > ";o2
9150 if o$ > "c" then print #3," " : print #3,"new power level (watts) = ";o2
9160 goto 9110
9170 if (o2 < 0 or o2 = 0) then o2 = pwr
9180 rem ----- ratio of power levels
9190 f1 = sqr(o2/pwr)
9200 if n$ = "h" then f1 = f1/s0/4/p
9210 print
9220 rem ----- designation of output file for near field data
9230 input "save to a file (y/n) > ";s$
9240 if s$ = "n" then 9320
9250 if s$ <> "y" then 9230
9260 input "filename (name.out) > ";f$
9270 if left$(right$(f$,4),1) = "." then 9280 : else f$ = f$+".out"
9280 if o$ > "c" then print #3," " : print #3,"filename (name.out) ";f$
9290 open f$ for output as #2
9300 print #2,nx*ny*nz;",";o2;",";n$
9310 rem ----- loop over z dimension
9320 for iz = 1 to nz
9330 zz = zi+(iz-1)*zc
9340 rem ----- loop over y dimension
9350 for iy = 1 to ny
9360 yy = yi+(iy-1)*yc
9370 rem ----- loop over x dimension
9380 for ix = 1 to nx
9390 xx = xi+(ix-1)*xc
9400 rem ----- near field header
9410 print #3," "
9420 if n$ = "e" then print #3,b$;"near electric fields";b$
9430 if n$ = "h" then print #3,b$;"near magnetic fields";b$
9440 print #3,tab (10);"field point: ";"x = ";xx;" y = ";yy;" z = ";zz
9450 print #3,"  vector","real","imaginary","magnitude","phase"
9460 if n$ = "e" then a$ = " v/m "
9470 if n$ = "h" then a$ = " amps/m "
9480 print #3," component  ",a$,a$,a$," deg"
9490 a1 = 0
9500 a3 = 0
9510 a4 = 0
9520 rem ----- loop over three vector components
9530 for i = 1 to 3
9540 x0 = xx
9550 y0 = yy
9560 z0 = zz
9570 if n$ = "h" then 9670
9580 t5 = 0
9590 t6 = 0
9600 t7 = 0
9610 if i = 1 then t5 = 2*s0
9620 if i = 2 then t6 = 2*s0
9630 if i = 3 then t7 = 2*s0
9640 u7 = 0
9650 u8 = 0
9660 goto 9770
9670 for j8 = 1 to 6
9680 k(j8,1) = 0
9690 k(j8,2) = 0
9700 next j8
9710 j9 = 1
9720 j8 = -1
9730 if i = 1 then x0 = xx+j8*s0/2
9740 if i = 2 then y0 = yy+j8*s0/2
9750 if i = 3 then z0 = zz+j8*s0/2
9760 rem ----- loop over source segments
9770 for j = 1 to n
9780 j1 = abs(c%(j,1))
9790 j2 = abs(c%(j,2))
9800 j3 = j2
9810 if j1 > j2 then j3 = j1
9820 f4 = sgn(c%(j,1))
9830 f5 = sgn(c%(j,2))
9840 f6 = 1
9850 f7 = 1
9860 u5 = 0
9870 u6 = 0
9880 rem ----- image loop
9890 for k = 1 to g step -2
9900 if c%(j,1) <> -c%(j,2) then 9960
9910 if k < 0 then 10570
9920 rem ----- compute vector potential a
9930 f6 = f4
9940 f7 = f5
9950 rem ----- compute psi(0,j,j+.5)
9960 p1 = 0
9970 p2 = 2*j3+j-1
9980 p3 = p2+0.5
9990 p4 = j2
10000 gosub 840
10010 u1 = t1*f5
10020 u2 = t2*f5
10030 rem ----- compute psi(0,j-.5,j)
10040 p3 = p2
10050 p2 = p2-0.5
10060 p4 = j1
10070 gosub 750
10080 v1 = f4*t1
10090 v2 = f4*t2
10100 rem ----- real part of vector potential contribution
10110 x3 = u1*ca(j2)+v1*ca(j1)
10120 y3 = u1*cb(j2)+v1*cb(j1)
10130 z3 = (f7*u1*cg(j2)+f6*v1*cg(j1))*k
10140 rem ----- imaginary part of vector potential contribution
10150 x5 = u2*ca(j2)+v2*ca(j1)
10160 y5 = u2*cb(j2)+v2*cb(j1)
10170 z5 = (f7*u2*cg(j2)+f6*v2*cg(j1))*k
10180 rem ----- magnetic field calculation completed
10190 if n$ = "h" then 10510
10200 d1 = (x3*t5+y3*t6+z3*t7)*w2
10210 d2 = (x5*t5+y5*t6+z5*t7)*w2
10220 rem ----- compute psi(.5,j,j+1)
10230 p1 = 0.5
10240 p2 = p3
10250 p3 = p3+1
10260 p4 = j2
10270 gosub 650
10280 u1 = t1
10290 u2 = t2
10300 rem ----- compute psi(-.5,j,j+1)
10310 p1 = -p1
10320 gosub 650
10330 u1 = (t1-u1)/s(j2)
10340 u2 = (t2-u2)/s(j2)
10350 rem ----- compute psi(.5,j-1,j)
10360 p1 = -p1
10370 p3 = p2
10380 p2 = p2-1
10390 p4 = j1
10400 gosub 650
10410 u3 = t1
10420 u4 = t2
10430 rem ----- compute psi(-.5,j-1,j)
10440 p1 = -p1
10450 gosub 650
10460 rem ----- gradient of scalar potential
10470 u5 = (u1+(u3-t1)/s(j1)+d1)*k+u5
10480 u6 = (u2+(u4-t2)/s(j1)+d2)*k+u6
10490 goto 10570
10500 rem ----- components of vector potential a
10510 k(1,j9) = k(1,j9)+(x3*cr(j)-x5*ci(j))*k
10520 k(2,j9) = k(2,j9)+(x5*cr(j)+x3*ci(j))*k
10530 k(3,j9) = k(3,j9)+(y3*cr(j)-y5*ci(j))*k
10540 k(4,j9) = k(4,j9)+(y5*cr(j)+y3*ci(j))*k
10550 k(5,j9) = k(5,j9)+(z3*cr(j)-z5*ci(j))*k
10560 k(6,j9) = k(6,j9)+(z5*cr(j)+z3*ci(j))*k
10570 next k
10580 if n$ = "h" then 10610
10590 u7 = u5*cr(j)-u6*ci(j)+u7
10600 u8 = u6*cr(j)+u5*ci(j)+u8
10610 next j
10620 if n$ = "e" then 10840
10630 rem ----- differences of vector potential a
10640 j8 = 1
10650 j9 = j9+1
10660 if j9 = 2 then 9730
10670 on i goto 10680,10730,10780
10680 h(3) = k(5,1)-k(5,2)
10690 h(4) = k(6,1)-k(6,2)
10700 h(5) = k(3,2)-k(3,1)
10710 h(6) = k(4,2)-k(4,1)
10720 goto 11060
10730 h(1) = k(5,2)-k(5,1)
10740 h(2) = k(6,2)-k(6,1)
10750 h(5) = h(5)-k(1,2)+k(1,1)
10760 h(6) = h(6)-k(2,2)+k(2,1)
10770 goto 11060
10780 h(1) = h(1)-k(3,2)+k(3,1)
10790 h(2) = h(2)-k(4,2)+k(4,1)
10800 h(3) = h(3)+k(1,2)-k(1,1)
10810 h(4) = h(4)+k(2,2)-k(2,1)
10820 goto 11060
10830 rem ----- imaginary part of electric field
10840 u7 = -m*u7/s0
10850 rem ----- real part of electric field
10860 u8 = m*u8/s0
10870 rem ----- magnitude and phase calculation
10880 s1 = 0
10890 if (u7 = 0 and u8 = 0) then 10910
10900 s1 = sqr(u7*u7+u8*u8)
10910 s2 = 0
10920 if u8 <> 0 then s2 = arctan(u7/u8)/p0
10930 if u8 > 0 then 10950
10940 s2 = s2+sgn(u7)*180
10950 if i = 1 then print #3,"   x  ",
10960 if i = 2 then print #3,"   y  ",
10970 if i = 3 then print #3,"   z  ",
10980 print #3,tab (15);f1*u8;tab (29);f1*u7;tab (43);f1*s1;tab (57);s2
10990 if s$ = "y" then print #2,f1*u8;",";f1*u7;",";f1*s1;",";s2
11000 rem ----- calculation for peak electric field
11010 s1 = s1*s1
11020 s2 = s2*p0
11030 a1 = a1+s1*cos(2*s2)
11040 a3 = a3+s1*sin(2*s2)
11050 a4 = a4+s1
11060 next i
11070 if n$ = "e" then 11300
11080 rem ----- magnetic field magnitude and phase calculation
11090 for i = 1 to 5 step 2
11100 s1 = 0
11110 if (h(i) = 0 and h(i+1) = 0) then 11130
11120 s1 = sqr(h(i)*h(i)+h(i+1)*h(i+1))
11130 s2 = 0
11140 if h(i) <> 0 then s2 = arctan(h(i+1)/h(i))/p0
11150 if h(i) > 0 then 11170
11160 s2 = s2+sgn(h(i+1))*180
11170 if i = 1 then print #3,"   x  ",
11180 if i = 3 then print #3,"   y  ",
11190 if i = 5 then print #3,"   z  ",
11200 print #3,tab (15);f1*h(i);tab (29);f1*h(i+1);tab (43);f1*s1;tab (57);s2
11210 if s$ = "y" then print #2,f1*h(i);",";f1*h(i+1);",";f1*s1;",";s2
11220 rem ----- calculation for peak magnetic field
11230 s1 = s1*s1
11240 s2 = s2*p0
11250 a1 = a1+s1*cos(2*s2)
11260 a3 = a3+s1*sin(2*s2)
11270 a4 = a4+s1
11280 next i
11290 rem ----- peak field calculation
11300 pk = sqr(a4/2+sqr(a1*a1+a3*a3)/2)
11310 print #3,"   maximum or peak field = ";f1*pk;a$
11320 if (s$ = "y" and n$ = "e") then print #2,f1*pk;",";o2
11330 if (s$ = "y" and n$ = "h") then print #2,f1*pk;",";o2
11340 if s$ = "y" then print #2,xx;",";yy;",";zz
11350 next ix
11360 next iy
11370 next iz
11380 close #2
11390 return
11400 rem ********** frequency input **********
11410 rem ----- set flag
11420 print
11430 input "frequency (MHz) > ";f
11440 if f = 0 then f = 299.8
11450 if o$ > "c" then print #3," " : print #3,"frequency (mhz):";f
11460 w = 299.8/f
11470 rem -----virtual dipole length for near field calculation
11480 s0 = 1.000000E-03*w
11490 rem ----- 1 / (4 * pi * omega * epsilon)
11500 m = 4.777834*w
11510 rem ----- set small radius modification condition
11520 srm = 1.000000E-04*w
11530 print #3,"    wave length = ";w;" meters"
11540 rem ----- 2 pi / wavelength
11550 w = 2*p/w
11560 w2 = w*w/2
11570 flg = 0
11580 return
11590 rem ********** geometry input **********
11600 rem ----- when geometry is changed, environment must be checked
11610 gosub 13780
11620 print
11630 if infile then 11690
11640 input "number of wires > ";nw
11650 if nw = 0 then return
11660 if nw <= mw then 11690
11670 print "number of wires exceeds dimension..."
11680 goto 11640
11690 if o$ > "c" then print #3," " : print #3,"no. of wires:";nw
11700 rem ----- initialize number of pulses to zero
11710 n = 0
11720 for i = 1 to nw
11730 if infile then gosub 15660 : goto 11990
11740 print
11750 print "wire no.";i
11760 input "   number of segments > ";s1
11770 if s1 = 0 then 11620
11780 a$ = "   end one coordinates (x,y,z) > "
11790 print a$;
11800 input x1,y1,z1
11810 if g < 0 and z1 < 0 then print "z cannot be negative" : goto 11790
11820 a$ = "   end two coordinates (x,y,z) > "
11830 print a$;
11840 input x2,y2,z2
11850 if g < 0 and z2 < 0 then print "z cannot be negative" : goto 11830
11860 if x1 = x2 and y1 = y2 and z1 = z2 then print "zero length wire." : goto 11750
11870 a$ = "   radius > "
11880 print "                     " a$;
11890 input a(i)
11900 if a(i) <= 0 then 11880
11910 rem ----- determine connections
11920 if o$ > "c" then print #3," " : print #3,"wire no.";i
11930 gosub 13080
11940 print "change wire no. ";i;" (y/n) > ";
11950 input a$
11960 if a$ = "y" then 11740
11970 if a$ <> "n" then 11940
11980 rem ----- compute direction cosines
11990 x3 = x2-x1
12000 y3 = y2-y1
12010 z3 = z2-z1
12020 d = sqr(x3*x3+y3*y3+z3*z3)
12030 ca(i) = x3/d
12040 cb(i) = y3/d
12050 cg(i) = z3/d
12060 s(i) = d/s1
12070 rem ----- compute connectivity data (pulses n1 to n)
12080 n1 = n+1
12090 n(i,1) = n1
12100 if (s1 = 1 and i1 = 0) then n(i,1) = 0
12110 n = n1+s1
12120 if i1 = 0 then n = n-1
12130 if i2 = 0 then n = n-1
12140 if n > mp then print "pulse number exceeds dimension" : close : goto 11640
12150 n(i,2) = n
12160 if (s1 = 1 and i2 = 0) then n(i,2) = 0
12170 if n < n1 then 12560
12180 for j = n1 to n
12190 c%(j,1) = i
12200 c%(j,2) = i
12210 w%(j) = i
12220 next j
12230 c%(n1,1) = i1
12240 c%(n,2) = i2
12250 rem ----- compute coordinates of break points
12260 i1 = n1+2*(i-1)
12270 i3 = i1
12280 x(i1) = x1
12290 y(i1) = y1
12300 z(i1) = z1
12310 if c%(n1,1) = 0 then 12390
12320 i2 = abs(c%(n1,1))
12330 f3 = sgn(c%(n1,1))*s(i2)
12340 x(i1) = x(i1)-f3*ca(i2)
12350 y(i1) = y(i1)-f3*cb(i2)
12360 if c%(n1,1) = -i then f3 = -f3
12370 z(i1) = z(i1)-f3*cg(i2)
12380 i3 = i3+1
12390 i6 = n+2*i
12400 for i4 = i1+1 to i6
12410 j = i4-i3
12420 x(i4) = x1+j*x3/s1
12430 y(i4) = y1+j*y3/s1
12440 z(i4) = z1+j*z3/s1
12450 next i4
12460 if c%(n,2) = 0 then 12540
12470 i2 = abs(c%(n,2))
12480 f3 = sgn(c%(n,2))*s(i2)
12490 i3 = i6-1
12500 x(i6) = x(i3)+f3*ca(i2)
12510 y(i6) = y(i3)+f3*cb(i2)
12520 if i = -c%(n,2) then f3 = -f3
12530 z(i6) = z(i3)+f3*cg(i2)
12540 goto 12640
12550 rem ---- single segmen 0 pulse case
12560 i1 = n1+2*(i-1)
12570 x(i1) = x1
12580 y(i1) = y1
12590 z(i1) = z1
12600 i1 = i1+1
12610 x(i1) = x2
12620 y(i1) = y2
12630 z(i1) = z2
12640 next i
12650 rem ********** geometry output **********
12660 print #3," "
12670 print #3,"                  **** antenna geometry ****"
12680 if n > 0 then 12730
12690 print
12700 print "number of pulses is zero....re-enter geometry"
12710 print
12720 goto 11640
12730 k = 1
12740 j = 0
12750 for i = 1 to n
12760 i1 = 2*w%(i)-1+i
12770 if k > nw then 12880
12780 if k = j then 12880
12790 j = k
12800 print #3," "
12810 print #3,"wire no. ";k;" coordinates",,,"connection pulse"
12820 print #3,"x","y","z","radius","end1 end2  no."
12830 if (n(k,1) <> 0 or n(k,2) <> 0) then 12880
12840 print #3,"-","-","-","    -"," -    -    0"
12850 k = k+1
12860 if k > nw then 12950
12870 goto 12790
12880 print #3,x(i1);tab (15);y(i1);tab (29);z(i1);tab (43);a(w%(i));tab (57);
12890 print #3,using "###  ###   ##";c%(i,1),c%(i,2),i
12900 if (i = n(k,2) or n(k,1) = n(k,2) or c%(i,2) = 0) then k = k+1
12910 if c%(i,1) = 0 then c%(i,1) = w%(i)
12920 if c%(i,2) = 0 then c%(i,2) = w%(i)
12930 if (k = nw and n(k,1) = 0 and n(k,2) = 0) then 12790
12940 if (i = n and k < nw) then 12790
12950 next i
12960 print
12970 close #1 : if infile then infile = 0 : if o$ > "c" then 13020
12980 input "    change geometry (y/n) > ";a$
12990 if a$ = "y" then 11620
13000 if a$ <> "n" then 12980
13010 rem ----- excitation input
13020 gosub 14390
13030 rem ----- loads/networks input
13040 gosub 14640
13050 flg = 0
13060 return
13070 rem ********** connections **********
13080 e(i) = x1
13090 l(i) = y1
13100 m(i) = z1
13110 e(i+nw) = x2
13120 l(i+nw) = y2
13130 m(i+nw) = z2
13140 g% = 0
13150 i1 = 0
13160 i2 = 0
13170 j1(i) = 0
13180 j2(i,1) = -i
13190 j2(i,2) = -i
13200 if g = 1 then 13320
13210 rem ----- check for ground connection
13220 if z1 = 0 then 13240
13230 goto 13270
13240 i1 = -i
13250 j1(i) = -1
13260 goto 13490
13270 if z2 = 0 then 13290
13280 goto 13320
13290 i2 = -i
13300 j1(i) = 1
13310 g% = 1
13320 if i = 1 then 13670
13330 for j = 1 to i-1
13340 rem ----- check for end1 to end1
13350 if (x1 = e(j) and y1 = l(j) and z1 = m(j)) then 13370
13360 goto 13420
13370 i1 = -j
13380 j2(i,1) = j
13390 if j2(j,1) = -j then j2(j,1) = j
13400 goto 13490
13410 rem ----- check for end1 to end2
13420 if (x1 = e(j+nw) and y1 = l(j+nw) and z1 = m(j+nw)) then 13440
13430 goto 13480
13440 i1 = j
13450 j2(i,1) = j
13460 if j2(j,2) = -j then j2(j,2) = j
13470 goto 13490
13480 next j
13490 if g% = 1 then 13670
13500 if i = 1 then 13670
13510 for j = 1 to i-1
13520 rem ----- check end2 to end2
13530 if (x2 = e(j+nw) and y2 = l(j+nw) and z2 = m(j+nw)) then 13550
13540 goto 13600
13550 i2 = -j
13560 j2(i,2) = j
13570 if j2(j,2) = -j then j2(j,2) = j
13580 goto 13670
13590 rem ----- check for end2 to end1
13600 if (x2 = e(j) and y2 = l(j) and z2 = m(j)) then 13620
13610 goto 13660
13620 i2 = j
13630 j2(i,2) = j
13640 if j2(j,1) = -j then j2(j,1) = j
13650 goto 13670
13660 next j
13670 print #3,"            coordinates","  ","  ","end         no. of"
13680 print #3,"   x","   y","   z","radius     connection     segments"
13690 print #3,x1;tab (15);y1;tab (29);z1;tab (57);i1
13700 print #3,x2;tab (15);y2;tab (29);z2;tab (43);a(i);tab (57);i2;tab (71);s1
13710 return
13720 rem ********** enviroment input **********
13730 print
13740 print "                        **** warning ****"
13750 print "redo geometry to ensure proper ground connection/disconnection"
13760 print
13770 rem ----- initialize number of radial wires to zero
13780 nr = 0
13790 rem ----- set environment
13800 print #3," "
13810 a$ = "environment (+1 for free space, -1 for ground plane) > "
13820 print a$;
13830 input g
13840 if o$ > "c" then print #3,a$;": ";g
13850 if g = 1 then 14370
13860 if g <> -1 then 13820
13870 rem ----- number of media
13880 a$ = " number of media (0 for perfectly conducting ground) > "
13890 print a$;
13900 input nm
13910 if nm <= mm then 13940
13920 print "number of media exceeds dimension..."
13930 goto 13890
13940 if o$ > "c" then print #3,a$;": ";nm
13950 rem ----- initialize boundary type
13960 tb = 1
13970 if nm = 0 then 14370
13980 if nm = 1 then 14050
13990 rem ----- type of boundary
14000 a$ = " type of boundary (1-linear, 2-circular) > "
14010 print "            ";a$;
14020 input tb
14030 if o$ > "c" then print #3,a$;": ";tb
14040 rem ----- boundary conditions
14050 for i = 1 to nm
14060 print "media";i
14070 a$ = " relative dielectric constant, conductivity > "
14080 print "         ";a$;
14090 input t(i),v(i)
14100 if o$ > "c" then print #3,a$;": ";t(i)"," v(i)
14110 if i > 1 then 14230
14120 if tb = 1 then 14230
14130 a$ = " number of radial wires in ground screen > "
14140 print "            ";a$;
14150 input nr
14160 if o$ > "c" then print #3,a$;": ";nr
14170 if nr = 0 then 14230
14180 a$ = " radius of radial wires > "
14190 print "                             ";a$;
14200 input rr
14210 if o$ > "c" then print #3,a$;": ";rr
14220 rem ----- initialize coordinate of media interface
14230 u(i) = 1000000
14240 rem ----- initialize height of media
14250 h(i) = 0
14260 if i = nm then 14310
14270 a$ = " x or r coordinate of next media interface > "
14280 print "          ";a$;
14290 input u(i)
14300 if o$ > "c" then print #3,a$;": ";u(i)
14310 if i = 1 then 14360
14320 a$ = " height of media > "
14330 print "                                    ";a$;
14340 input h(i)
14350 if o$ > "c" then print #3,a$;": ";h(i)
14360 next i
14370 return
14380 rem ********** excitation input **********
14390 print
14400 a$ = "no. of sources > "
14410 print a$;
14420 input ns
14430 if ns < 1 then ns = 1
14440 if ns <= mp then 14470
14450 print "no. of sources exceeds dimension ..."
14460 goto 14410
14470 if o$ > "c" then print #3," " : print #3,a$;": ";ns
14480 for i = 1 to ns
14490 print
14500 print "source no. ";i;":"
14510 a$ = "pulse no., voltage magnitude, phase (degrees) > "
14520 print a$;
14530 input e(i),vm,vp
14540 if e(i) <= n then 14570
14550 print "pulse number exceeds number of pulses..."
14560 goto 14520
14570 if o$ > "c" then print #3,a$;": ";e(i)"," vm "," vp
14580 l(i) = vm*cos(vp*p0)
14590 m(i) = vm*sin(vp*p0)
14600 next i
14610 if flg = 2 then flg = 1
14620 return
14630 rem ********** loads input **********
14640 print
14650 input "number of loads      > ";nl
14660 if nl <= ml then 14690
14670 print "number of loads exceeds dimension..."
14680 goto 14650
14690 if o$ > "c" then print #3,"number of loads";nl
14700 if nl < 1 then 15010
14710 input "s-parameter (s=jw) impedance load (y/n) > ";l$
14720 if l$ <> "y" and l$ <> "n" then 14710
14730 a$ = "pulse no.,resistance,reactance"
14740 if l$ = "y" then a$ = "pulse no., order of s-parameter function"
14750 for i = 1 to nl
14760 print
14770 print "load no. ";i;": "
14780 if l$ = "y" then 14850
14790 print a$;
14800 input lp(i),la(1,i,1),la(2,i,1)
14810 if lp(i) > n then print "pulse number exceeds number of pulses..." : goto 14790
14820 if o$ > "c" then print #3,a$;": ";lp(i);",";la(1,i,1);",";la(2,i,1)
14830 goto 15000
14840 rem ----- s-parameter loads
14850 print a$;" > ";
14860 input lp(i),ls(i)
14870 if lp(i) > n then print "pulse number exceeds number of pulses..." : goto 14850
14880 if ls(i) > ma then print "maximum dimension is 10" : goto 14860
14890 if o$ > "c" then print #3,a$;": ";lp(i);",";ls(i)
14900 for j = 0 to ls(i)
14910 a$ = "numerator, denominator coefficients of s^"
14920 print a$;j;" > ";
14930 input la(1,i,j),la(2,i,j)
14940 if o$ > "c" then print #3,a$;j;":";la(1,i,j);",";la(2,i,j)
14950 next j
14960 if ls(i) > 0 then 15000
14970 ls(i) = 1
14980 la(1,i,1) = 0
14990 la(2,i,1) = 0
15000 next i
15010 flg = 0
15020 return
15030 rem ********** main program **********
15040 rem ----- data initialization
15050 rem ----- pi
15060 p = 4*arctan(1)
15070 rem ----- changes degrees to radians
15080 p0 = p/180
15090 b$ = "********************"
15100 rem ----- intrinsic impedance of free space divided by 2 pi
15110 g0 = 29.979221
15120 rem ---------- q-vector for gaussian quadrature
15130 read q(1),q(2),q(3),q(4),q(5),q(6),q(7),q(8),q(9),q(10),q(11),q(12)
15140 read q(13),q(14)
15150 data 0.288675,0.5,0.430568,0.173927,0.169991,0.326073
15160 data 0.480145,0.050614,0.398333,0.111191
15170 data 0.262766,0.156853,0.091717,0.181342
15180 rem ---------- e-vector for coefficients of elliptic integral
15190 read c0,c1,c2,c3,c4,c5,c6,c7,c8,c9
15200 data 1.386294,0.096663,0.035901,0.037426,0.014512
15210 data 0.5,0.124986,0.068802,0.033284,4.417870E-03
15220 rem ----- identify output device
15230 gosub 15890
15240 print #3,tab (20);b$;b$
15250 print #3,tab (22);"mini-numerical electromagnetics code"
15260 print #3,tab (36);"mininec"
15270 print #3,tab (24);date$;tab (48);time$
15280 print #3,tab (20);b$;b$
15290 rem ----- frequency input
15300 gosub 11420
15310 rem ----- environment input
15320 gosub 13780
15330 rem ----- check for nec-type geometry input
15340 gosub 15590
15350 rem ----- geometry input
15360 gosub 11620
15370 rem ----- menu
15380 print
15390 print b$;"    mininec menu    ";b$
15400 print "   g - change geometry     c - compute/display currents"
15410 print "   e - change environment  p - compute far-field patterns"
15420 print "   x - change excitation   n - compute near-fields"
15430 print "   l - change loads"
15440 print "   f - change frequency    q - quit"
15450 print b$;b$;b$
15460 input "   command > ";c$
15470 if c$ = "f" then gosub 11420
15480 if c$ = "p" then gosub 6300
15490 if c$ = "x" then gosub 14390
15500 if c$ = "e" then gosub 13730
15510 if c$ = "g" then gosub 11610
15520 if c$ = "c" then gosub 5060
15530 if c$ = "l" then gosub 14640
15540 if c$ = "n" then gosub 8840
15541 if len(c$) = 0 then bye
15550 if c$ <> "q" then 15380
15560 if o$ = "p" then print #3,chr$(12) : else if o$ = "c" then print #3," "
15570 rem close 
15571 close #1 : close #2 : close #3 : close #4
15580 goto 16260
15590 rem ********** nec-type geometry input **********
15600 rem open "mininec.inp" as #1 len = 30
15610 rem field #1,2 as s$,4 as x1$,4 as y1$,4 as z1$,4 as x2$,4 as y2$,4 as z2$,4 as r$
15612 nw = 0 : infile = 0
15613 open "mininec.inp" for input as #1 else goto 15640
15614 print "mininec.inp opened"
15620 rem get 1
15621 input #1, tmp1$
15630 rem nw = cvi(s$)
15631 nw = val(field$(tmp1$,1))
15640 if nw then infile = 1
15641 print "infile = ",infile
15650 return
15660 rem ---------- get geometry data from mininec.inp
15670 rem get 1
15672 input #1, tmp1$
15673 print "mininec.inp line = ", tmp1$
15680 rem s1 = cvi(s$)
15681 s1 = val(field$(tmp1$,1))
15690 rem x1 = cvs(x1$)
15691 x1 = val(field$(tmp1$,2))
15700 rem y1 = cvs(y1$)
15701 y1 = val(field$(tmp1$,3))
15710 rem z1 = cvs(z1$)
15711 z1 = val(field$(tmp1$,4))
15720 rem x2 = cvs(x2$)
15721 x2 = val(field$(tmp1$,5))
15730 rem y2 = cvs(y2$)
15731 y2 = val(field$(tmp1$,6))
15740 rem z2 = cvs(z2$)
15741 z2 = val(field$(tmp1$,7))
15750 rem a(i) = cvs(r$)
15751 a(i) = val(field$(tmp1$,8))
15752 print s1,x1,y1,z1,x2,y2,z2,a(i)
15760 if g < 0 then if z1 < 0 or z2 < 0 then gosub 15810
15770 print #3," " : print #3,"wire no.";i
15780 if x1 = x2 and y1 = y2 and z1 = z2 then print "wire length is zero." : goto 15560
15790 gosub 13080
15800 return
15810 if izneg then 15850
15820 print "negative z value encountered for ground plane."
15830 input "abort or convert negative z value to zero (a/c)? > ";a$
15840 rem if a$ = "a" then 15560 : else if a$ = "c" then izneg = 1 : else 15830
15841 if a$ = "a" then goto 15560 
15842 if a$ = "c" then izneg = 1 : else goto 15830
15850 if z1 < 0 then z1 = -z1
15860 if z2 < 0 then z2 = -z2
15870 return
15880 rem ********** identify output device **********
15890 input "output to console, printer, or disk (c/p/d) > ";o$
15900 if o$ = "c" then f$ = "stdout" : goto 15950
15910 if o$ = "p" then f$ = "lpt1:" : goto 15950
15920 if o$ <> "d" then 15890
15930 input "filename (name.out) > ";f$
15940 if left$(right$(f$,4),1) = "." then 15950 : else f$ = f$+".out"
15950 open f$ for output as #3
15960 cls
15970 return
15980 rem ********** calculate elapsed time **********
15990 ih = val(mid$(t$,1,2))-val(mid$(ot$,1,2))
16000 im = val(mid$(t$,4,2))-val(mid$(ot$,4,2))
16010 is0 = val(mid$(t$,7,2))-val(mid$(ot$,7,2))
16020 if is0 < 0 then is0 = is0 +60 : im = im-1
16030 if im < 0 then im = im+60 : ih = ih-1
16040 if ih < 0 then ih = ih+24
16050 t$ = ":"+mid$(str$(is0 +100),3)
16060 if ih then t$ = mid$(str$(ih),2)+":"+mid$(str$(im+100),3)+t$ : else t$ = mid$(str$(im),2)+t$
16070 return
16080 rem ********** calculate approximate time remaining **********
16090 ipct = 100*pct
16100 t$ = time$
16110 ih = val(mid$(t$,1,2))-val(mid$(ot$,1,2))
16120 if ih < 0 then ih = ih+24
16130 im = val(mid$(t$,4,2))-val(mid$(ot$,4,2))
16140 is0 = val(mid$(t$,7,2))-val(mid$(ot$,7,2))
16150 is0 = is0 +60*(im+60*ih)
16160 is0 = is0 *(1/pct-1)
16170 im = int(is0 /60)
16180 is0 = is0 mod 60
16190 ih = int(im/60)
16200 im = im mod 60
16210 t$ = ":"+mid$(str$(is0 +100),3)
16220 if ih then t$ = mid$(str$(ih),2)+":"+mid$(str$(im+100),3)+t$ : else t$ = mid$(str$(im),2)+t$
16230 rem locate csrlin,1
16231 cls : gotoxy csrlin,1
16240 print q$;ipct;"% complete - approx time remaining " t$ "   ";
16250 return
16260 bye
16270 end
99990 end
