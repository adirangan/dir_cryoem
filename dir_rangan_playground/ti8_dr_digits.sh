nMode=0
while [ $nMode -lt 12 ]
do
   ./ti8_dr_digits.exe  $nMode
   if [ $nMode -eq 12 ]
   then
      break
   fi
   nMode=`expr $nMode + 1`
done

