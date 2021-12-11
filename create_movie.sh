#!/bin/sh

Nframes=2625
fname="frame_data.dat"
povfile="magnet"
#gifanim=/home/isouki/apps/whirlgif/whirlgif
#animtag="-o movie.gif "

j=0
while [ "${j}" -le "${Nframes}" ]
do
mv ${j}.dat ${fname}
povray -D +A +H1600 -J +Q9 +R5 -V +W1600 ${povfile}.pov
mv ${povfile}.png ${j}.png
mv ${fname} ${j}.dat
#convert -resize 800x800 ${j}.png small_${j}.png
j=$(($j+1))
done

ffmpeg -i '%d.png' -s 800x800 -r 30 -vcodec qtrle movie.mov
#ffmpeg -f image2 -i small_%d.png movie.mpg

#for (( i=1;i<=${Nframes};i++ ))
#do
#    animtag="${animtag} -time 3 ${i}.gif "
#done
#whirlgif ${animtag}

