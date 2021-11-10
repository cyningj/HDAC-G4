for i in $(cat list)
do
        {
                python screen.py $i -e -w 60 -s 5 > $i"G4-clips"
}&
done
wait
echo 'success!'
wait
head -n 1 list | while read line
do
        mkdir ./$line"done"/
        mv *_G4 ./$line"done"/
        mv list ./$line"done"/
done
