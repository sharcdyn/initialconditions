program xyz2geom
implicit none
integer nat
character*2, allocatable :: at(:)
real*8, allocatable :: geom(:,:)
real*8, allocatable :: noat(:),mass(:)
character*100 kkchar

real*8 :: bohr2ang=0.52917721067

integer i,j

write(6,*) "XYZ file to be converted"
read(5,*) kkchar
open(1,file=kkchar,status="old",iostat=i)
if (i.ne.0) stop "File does not exist"
read(1,*) nat
read(1,*)
allocate(at(nat),geom(nat,3),noat(nat),mass(nat))
do i=1,nat
 read(1,*) at(i),(geom(i,j),j=1,3)
 do j=1,3
  geom(i,j)=geom(i,j)/bohr2ang
 enddo
 call getnoatandmass(at(i),noat(i),mass(i))
enddo

write(6,*) "Select type of output file"
read(5,*) kkchar
select case(kkchar)
 case("geom")
  write(6,*) "Creating file geom"
  open(1,file="geom")
  do i=1,nat
   write(1,"(A2,x,F5.1,4(x,F30.20))") at(i),noat(i),geom(i,:),mass(i)
  enddo
 case("bagel")
  write(6,*) "Creating file geom.json"
  open(1,file="geom.json")
  write(1,"(A)") '{ "title" : "molecule", '
  write(1,"(A)") '"basis" : "",'
  write(1,"(A)") '"df_basis" : "",'
  write(1,"(A)") '"angstrom" : true,'
  write(1,"(A)") '"geometry" : [ '
  do i=1,nat
   if (i.ne.nat) then
    write(1,"(A,A,A,4(x,F20.10,A))") &
     '{"atom" : "',trim(at(i)),'","xyz" : [ ',&
      geom(i,1),', ',geom(i,2),', ',geom(i,3),&
     ' ], "mass" : ',mass(i),' },'
   else
    write(1,"(A,A,A,4(x,F20.10,A))") &
     '{"atom" : "',trim(at(i)),'","xyz" : [ ',&
      geom(i,1),', ',geom(i,2),', ',geom(i,3),&
     ' ], "mass" : ',mass(i),' }] }'
   endif
  enddo
 case default
  write(6,*) "Valid output files are geom or bagel"
end select

end

subroutine getnoatandmass(at,noat,mass)
implicit none
character*2, intent(in) :: at
real*8, intent(out) :: noat,mass

integer i
character*2 sat(118)
real*8 atmass(118)

sat(1)="H"
sat(2)="He"
sat(3)="Li"
sat(4)="Be"
sat(5)="B"
sat(6)="C"
sat(7)="N"
sat(8)="O"
sat(9)="F"
sat(10)="Ne"
sat(11)="Na"
sat(12)="Mg"
sat(13)="Al"
sat(14)="Si"
sat(15)="P"
sat(16)="S"
sat(17)="Cl"
sat(18)="Ar"
sat(19)="K"
sat(20)="Ca"
sat(21)="Sc"
sat(22)="Ti"
sat(23)="V"
sat(24)="Cr"
sat(25)="Mn"
sat(26)="Fe"
sat(27)="Co"
sat(28)="Ni"
sat(29)="Cu"
sat(30)="Zn"
sat(31)="Ga"
sat(32)="Ge"
sat(33)="As"
sat(34)="Se"
sat(35)="Br"
sat(36)="Kr"
sat(37)="Rb"
sat(38)="Sr"
sat(39)="Y"
sat(40)="Zr"
sat(41)="Nb"
sat(42)="Mo"
sat(43)="Tc"
sat(44)="Ru"
sat(45)="Rh"
sat(46)="Pd"
sat(47)="Ag"
sat(48)="Cd"
sat(49)="In"
sat(50)="Sn"
sat(51)="Sb"
sat(52)="Te"
sat(53)="I"
sat(54)="Xe"
sat(55)="Cs"
sat(56)="Ba"
sat(57)="La"
sat(58)="Ce"
sat(59)="Pr"
sat(60)="Nd"
sat(61)="Pm"
sat(62)="Sm"
sat(63)="Eu"
sat(64)="Gd"
sat(65)="Tb"
sat(66)="Dy"
sat(67)="Ho"
sat(68)="Er"
sat(69)="Tm"
sat(70)="Yb"
sat(71)="Lu"
sat(72)="Hf"
sat(73)="Ta"
sat(74)="W"
sat(75)="Re"
sat(76)="Os"
sat(77)="Ir"
sat(78)="Pt"
sat(79)="Au"
sat(80)="Hg"
sat(81)="Tl"
sat(82)="Pb"
sat(83)="Bi"
sat(84)="Po"
sat(85)="At"
sat(86)="Rn"
sat(87)="Fr"
sat(88)="Ra"
sat(89)="Ac"
sat(90)="Th"
sat(91)="Pa"
sat(92)="U"
sat(93)="Np"
sat(94)="Pu"
sat(95)="Am"
sat(96)="Cm"
sat(97)="Bk"
sat(98)="Cf"
sat(99)="Es"
sat(100)="Fm"
sat(101)="Md"
sat(102)="No"
sat(103)="Lr"
sat(104)="Rf"
sat(105)="Db"
sat(106)="Sg"
sat(107)="Bh"
sat(108)="Hs"
sat(109)="Mt"
sat(110)="Ds"
sat(111)="Rg"
sat(112)="Cn"
sat(113)="Nh"
sat(114)="Fl"
sat(115)="Mc"
sat(116)="Lv"
sat(117)="Ts"
sat(118)="Og"
atmass(1)=1.00782503223
atmass(2)=3.0160293201
atmass(3)=6.0151228874
atmass(4)=9.012183065
atmass(5)=10.01293695
atmass(6)=12.0000000
atmass(7)=14.00307400443
atmass(8)=15.99491461957
atmass(9)=18.99840316273
atmass(10)=19.9924401762
atmass(11)=22.9897692820
atmass(12)=23.985041697
atmass(13)=26.98153853
atmass(14)=27.97692653465
atmass(15)=30.97376199842
atmass(16)=31.9720711744
atmass(17)=34.968852682
atmass(18)=35.967545105
atmass(19)=38.9637064864
atmass(20)=39.962590863
atmass(21)=44.95590828
atmass(22)=45.95262772
atmass(23)=49.94715601
atmass(24)=49.94604183
atmass(25)=54.93804391
atmass(26)=53.93960899
atmass(27)=58.93319429
atmass(28)=57.93534241
atmass(29)=62.92959772
atmass(30)=63.92914201
atmass(31)=68.9255735
atmass(32)=69.92424875
atmass(33)=74.92159457
atmass(34)=73.922475934
atmass(35)=78.9183376
atmass(36)=77.92036494
atmass(37)=84.9117897379
atmass(38)=83.9134191
atmass(39)=88.9058403
atmass(40)=89.9046977
atmass(41)=92.9063730
atmass(42)=91.90680796
atmass(43)=96.9063667
atmass(44)=95.90759025
atmass(45)=102.9054980
atmass(46)=101.9056022
atmass(47)=106.9050916
atmass(48)=105.9064599
atmass(49)=112.90406184
atmass(50)=111.90482387
atmass(51)=120.9038120
atmass(52)=119.9040593
atmass(53)=126.9044719
atmass(54)=123.9058920
atmass(55)=132.9054519610
atmass(56)=129.9063207
atmass(57)=137.9071149
atmass(58)=135.90712921
atmass(59)=140.9076576
atmass(60)=141.9077290
atmass(61)=144.9127559
atmass(62)=143.9120065
atmass(63)=150.9198578
atmass(64)=151.9197995
atmass(65)=158.9253547
atmass(66)=155.9242847
atmass(67)=164.9303288
atmass(68)=161.9287884
atmass(69)=168.9342179
atmass(70)=167.9338896
atmass(71)=174.9407752
atmass(72)=173.9400461
atmass(73)=179.9474648
atmass(74)=179.9467108
atmass(75)=184.9529545
atmass(76)=183.9524885
atmass(77)=190.9605893
atmass(78)=189.9599297
atmass(79)=196.96656879
atmass(80)=195.9658326
atmass(81)=202.9723446
atmass(82)=203.9730440
atmass(83)=208.9803991
atmass(84)=208.9824308
atmass(85)=209.9871479
atmass(86)=210.9906011
atmass(87)=223.0197360
atmass(88)=223.0185023
atmass(89)=227.0277523
atmass(90)=230.0331341
atmass(91)=231.0358842
atmass(92)=233.0396355
atmass(93)=236.046570
atmass(94)=238.0495601
atmass(95)=241.0568293
atmass(96)=243.0613893
atmass(97)=247.0703073
atmass(98)=249.0748539
atmass(99)=252.082980
atmass(100)=257.0951061
atmass(101)=258.0984315
atmass(102)=259.10103
atmass(103)=262.10961
atmass(104)=267.12179
atmass(105)=268.12567
atmass(106)=271.13393
atmass(107)=272.13826
atmass(108)=270.13429
atmass(109)=276.15159
atmass(110)=281.16451
atmass(111)=280.16514
atmass(112)=285.17712
atmass(113)=284.17873
atmass(114)=289.19042
atmass(115)=288.19274
atmass(116)=293.20449
atmass(117)=292.20746
atmass(118)=294.21392

do i=1,118
 if (at==sat(i)) then
  noat=float(i)
  mass=atmass(i)
 endif
enddo


end

