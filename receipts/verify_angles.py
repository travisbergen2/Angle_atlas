"""Recompute every derived angle/length quoted in the angle-and-length list (2026-09-09).
Each line prints the value; asserts pin the values quoted in the thread records."""
import math
d = math.degrees
pi = math.pi
g1 = 14.134725141734693  # gamma_1 (Odlyzko)
g100k = 74920.827498994  # gamma_100000 (Odlyzko table top, per 09-04 record)

print("== 1. the dial theta(n) = 90/2^Omega ==")
for k in range(6): print(f"Omega={k}: {90/2**k}°")
print("== 2. exponent-lattice angle of 2^k p from the 2-axis: arctan(1/k) ==")
lat = [90.0] + [d(math.atan(1/k)) for k in range(1,7)]
print([round(x,3) for x in lat])
assert abs(lat[2]-26.565)<1e-3 and abs(lat[3]-18.435)<1e-3 and abs(lat[4]-14.036)<1e-3
print("== 3. shelf heights c(m) = m(pi/2)/2^Omega(m) ==")
def Omega(n):
    c=0; p=2
    while p*p<=n:
        while n%p==0: n//=p; c+=1
        p+=1
    return c+(1 if n>1 else 0)
for m in (1,3,5,7,9,15): print(m, round(m*(pi/2)/2**Omega(m),4))
assert abs(9*(pi/2)/2**Omega(9)-3.5343)<1e-4
print("== 4. nested radicals 2cos(90/2^k), Viete ==")
print([round(2*math.cos(math.radians(90/2**k)),5) for k in range(5)])
prod=1.0
for k in range(1,60): prod*=math.cos(math.radians(90/2**k))
print("Viete product", round(prod,10), "2/pi", round(2/pi,10)); assert abs(prod-2/pi)<1e-9
print("== 5. Gaussian angles ==")
print("arg(2+i) =", round(d(math.atan(1/2)),3))
print("angle between (2,1),(2,-1) = arccos(3/5) =", round(d(math.acos(3/5)),2))
print("arg(3+4i) = arctan(4/3) =", round(d(math.atan(4/3)),2))
print("(1,1).(1,-1) =", 1*1+1*(-1), "| lengths", round(math.sqrt(2),4), round(math.sqrt(2),4), "sum length", abs(complex(1,1)+complex(1,-1)))
print("== 6. (2,3,4) triangle largest angle ==")
C = d(math.acos((4+9-16)/(2*2*3))); print(round(C,2)); assert abs(C-104.48)<0.01
print("== 7. arctan((p+1)/p) drift ==")
for p in (2,3,5,7,11,101): print(p, round(d(math.atan((p+1)/p)),2))
assert abs(d(math.atan(3/2))-56.31)<0.01 and abs(d(math.atan(4/3))-53.13)<0.01 and abs(d(math.atan(6/5))-50.19)<0.01
print("== 8. balance angle theta = arctan(beta/(1-beta)) ==")
for b in (0.35,0.5): 
    t=d(math.atan(b/(1-b))); print(b, round(t,2), round(90-t,2))
assert abs(d(math.atan(0.35/0.65))-28.30)<0.01
print("== 9. sopfr dial arctan(sopfr(n)/n) ==")
def sopfr(n):
    s=0; p=2
    while p*p<=n:
        while n%p==0: n//=p; s+=p
        p+=1
    return s+(n if n>1 else 0)
for n in (3,4,6,45,128,3087): print(n, sopfr(n), round(d(math.atan(sopfr(n)/n)),3))
assert abs(d(math.atan(5/6))-39.81)<0.01 and abs(d(math.atan(11/45))-13.736)<0.001 and abs(d(math.atan(27/3087))-0.501)<0.001
print("alpha^-1/10 =", round(137.035999/10,3))
print("== 10. Bergen Disc (Li map z = 1 - 1/s, R = 400 px) ==")
def rim(g): return d(math.atan2(g, g*g-0.25))
print("rim angle gamma_1:", round(rim(g1),3), "| 1/gamma_1 in deg:", round(d(1/g1),3))
print("rim angle #100000:", round(rim(g100k),6))
print("rim span per side (100k zeros):", round(rim(g1)-rim(g100k),2)); assert abs(rim(g1)-4.05)<0.01
print("arg(rho_1) =", round(d(math.atan2(g1,0.5)),2)); assert abs(d(math.atan2(g1,0.5))-87.97)<0.01
print("arg(rho_100000) =", round(d(math.atan2(g100k,0.5)),4))
print("zero #100000 distance to vanishing point (px):", round(400/math.hypot(0.5,g100k),4))
print("off-rim offset of delta=0.1 at gamma_1 (px):", round(400*0.1/(0.6**2+g1**2),2), "| sub-pixel for gamma >", round(math.sqrt(40),2))
print("Weil-leg flip threshold delta* = pi/(4 a^2 gamma_1): a=1:", round(pi/(4*g1),4), " a=0.5:", round(pi/(4*0.25*g1),4))
print("phantom 0.7+30i bearing from 1/2:", round(d(math.atan2(30,0.2)),3))
print("== 11. anti-numerology exhibits ==")
print("10*sqrt2 =", round(10*math.sqrt(2),4), "gamma_1 =", round(g1,4), "miss", round(10*math.sqrt(2)-g1,4))
print("== 12. Collatz dial ==")
print("Omega(27) =", Omega(27), "angle", 90/2**Omega(27), "r*theta =", 27*90/2**Omega(27))
print("log2(3/4) =", round(math.log2(3/4),5))
print("== 13. sector windows (E-RACE-A) in degrees ==")
print({k: round(d(v),1) for k,v in {"pi/8":pi/8,"pi/4":pi/4,"pi/2":pi/2,"pi":pi,"3pi/4":3*pi/4}.items()})
print("== 14. Lorentz factor ==")
for v in (0.6,0.8): g=1/math.sqrt(1-v*v); print(f"v={v}c: gamma={g:.3f}; moving rod L=L0/{g:.3f}; moving clock dt={g:.3f} dtau")
print("ALL ASSERTS PASSED")
