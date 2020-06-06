#!/bin/bash -l

skip_col=2
if [ $1 == '-col' ]; then
   skip_col=$2
   shift 2
fi


head -n 10 $1
echo "#"
echo "# Results based on $* files"
echo "#"
echo "#LGCMC output column sequence:"
sed -n '11,12p' $1
echo "#This script will not calculate mean and sd for first $skip_col colums."
echo "#After that first the mean will be printed following LGCMC column sequence."
echo "#Then relative standard deviation (in %) will be printed following LGCMC column sequence."
echo "#"

awk -v skipcol=$skip_col '{ if (FNR==1) file++
  for (i=1;i<=NF;i++) {
    sum[FNR" "i]+=$i
    count[FNR" "i]++
    data[file" "FNR" "i]=$i
  }
}END{
  skipline=12                #No of line to skip
                  
  for (i=skipline+1;i<=FNR;i++) {
    for (j=1;j<=NF;j++) {
      avg[i" "j]=sum[i" "j]/count[i" "j]
    }
  }
  for (f=1;f<=file;f++) {
    for (i=skipline+1;i<=FNR;i++) {
      for (j=1;j<=NF;j++) {
        tmp[i" "j]+=(data[f" "i" "j]-avg[i" "j])*(data[f" "i" "j]-avg[i" "j]);
      }
    }
  }
  for (i=skipline+1;i<=FNR;i++) {
    for (j=1;j<=NF;j++) {
      printf "%18.8g", avg[i" "j]

    }
    for (j=skipcol+1;j<=NF;j++) {
#      printf "%18.8g", sqrt(tmp[i" "j]/count[i" "j]) 
      if(avg[i" "j]!=0)
      printf "%18.8g", sqrt(tmp[i" "j]/count[i" "j])/avg[i" "j]*100 
      else
      printf "         NAN      "
    }

    printf "\n"
  }
}' $* 


