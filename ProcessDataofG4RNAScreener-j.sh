#!/bin/bash
#去掉标题并将chrx:xxx-XXX转化为chr      xxx-XXX转化为chr        xxx     xxx
ls > list
sed -i 's/list//g'  list
sed -i 's/ProcessDataofG4RNAScreener-j.sh//g'  list
for i in `cat list`
do
        awk -F"\t" 'BEGIN{OFS="\t"}{gsub(/-/,"\t",$2);gsub(/:/,"\t",$2);print $0}' $i > $i"_clean"
done
