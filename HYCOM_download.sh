! This is a bash and ferret script that downloads and saves daily hycom files 

! We are setting a region of 120W:120E 0:50N and starting this on Jan 1, 2000, which is time stamp 2648
! The files need to have the modulo changed to 360 before saving and we are calculating an averaged depth of 1:25, 0:50, 0:100 in addition to having the top 50 m water layers


use http://apdrc.soest.hawaii.edu/erddap/griddap/hawaii_soest_3151_b22a_aae6
set axis/modulo=360 longitude
set region/x=120:240/y=0:50

let avg_U_25 = water_u[k=1:10@ave]
let avg_U_50 = water_u[k=1:15@ave]
let avg_U_100 = water_u[k=1:20@ave]
let avg_V_25 = water_v[k=1:10@ave]
let avg_V_50 = water_v[k=1:15@ave]
let avg_V_100 = water_v[k=1:20@ave]

repeat/L=2648:2652 (save/clobber/file=HYCOM_LAPS_`l`.nc/L=`l` water_u[L=`l`, k=1], avg_U_25[L=`l`], avg_U_50[L=`l`], avg_U_100[L=`l`], water_u[L=`l`, k=1], avg_V_25[L=`l`], avg_V_50[L=`l`], avg_V_100[L=`l`])
