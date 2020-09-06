#!/bin/bash
cd '/..'
R CMD BATCH --vanilla '--args m=1 f=1 adv_id=1 age.truncate.female=18 l_m=7 K_m=5 l_y=11 K_y=7' main_template.R main_1.out &
R CMD BATCH --vanilla '--args m=1 f=1 adv_id=2 age.truncate.female=18 l_m=7 K_m=5 l_y=11 K_y=7' main_template.R main_2.out &
R CMD BATCH --vanilla '--args m=1 f=1 adv_id=3 age.truncate.female=18 l_m=7 K_m=5 l_y=11 K_y=7' main_template.R main_3.out &
R CMD BATCH --vanilla '--args m=1 f=1 adv_id=4 age.truncate.female=18 l_m=7 K_m=5 l_y=11 K_y=7' main_template.R main_4.out &
R CMD BATCH --vanilla '--args m=1 f=1 adv_id=5 age.truncate.female=18 l_m=7 K_m=5 l_y=11 K_y=7' main_template.R main_5.out &
R CMD BATCH --vanilla '--args m=1 f=1 adv_id=6 age.truncate.female=18 l_m=7 K_m=5 l_y=11 K_y=7' main_template.R main_6.out &
R CMD BATCH --vanilla '--args m=1 f=1 adv_id=8 age.truncate.female=18 l_m=7 K_m=5 l_y=11 K_y=7' main_template.R main_8.out &
R CMD BATCH --vanilla '--args m=1 f=1 adv_id=9 age.truncate.female=18 l_m=7 K_m=5 l_y=11 K_y=7' main_template.R main_9.out &
