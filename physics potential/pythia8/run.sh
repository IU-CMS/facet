#!bin/bash -e

#----------
#if on, initialize; else, don't.
pythia=on
lhe2root=on
p4_filter=on
#----------

dir="$PWD"
sub_dir="test2"
if [[ ! -d "$sub_dir" ]]; then mkdir ${sub_dir}; echo "$sub_dir"" is created."; fi

pythia_dir=~/"Programs/pythia8244/examples"
ExRoot_dir=~/"Programs/MG5_aMC_v2_6_5/ExRootAnalysis"
unique_name="HeavyMesonProduction"
#ExRoot_dir="/opt/MG5_aMC_v2_7_2/ExRootAnalysis/"
ncore=4
nEvents=1000 # each thread
eCM=13000
# particle names in filter.C should be edited.
pid=443
par="J_Psy"
four_vec="pyt8_10MEvt_13TeVpp_${par}_4Vec_${sub_dir}"

#--------------------------------------------

pythia(){

cd $pythia_dir

lhe_file="$dir"/"$sub_dir"/"$unique_name"
log_file="$lhe_file"


pTHatMin=(2)
pTHatMax=(50)

if [[ ${#pTHatMin[@]} -ne ${#pTHatMax[@]} ]]; then
        echo "mission aborted"
        exit
fi

for (( i = 0; i < ncore ; i++ ));
do

time ./upsilonProduction "$lhe_file"_"$i".lhe $eCM ${pTHatMin[0]} ${pTHatMax[0]} $RANDOM $nEvents > "$log_file"_"$i".dat &

done
wait

cd $dir
}



lhe2root(){

cd $ExRoot_dir

for(( i=0 ; i<ncore ; i++ ))
do
	if [ -f "$dir"/${sub_dir}/"$unique_name"_"$i".lhe ]; then
		./ExRootLHEFConverter "$dir"/"${sub_dir}"/"$unique_name"_"$i".lhe "$dir"/"${sub_dir}"/"$unique_name"_"$i".root &
	fi
done
wait

cd $dir
}



p4_filter(){

root_list=()
root_four_vec=()


for (( i = 0; i < ncore; i++ ))
do

	if [ -f "${sub_dir}"/"${unique_name}"_"$i".root ]; then

		root_list+=("$unique_name"_"$i".root)
	
		root_four_vec+=(${sub_dir}/${four_vec}_${i}.root)
	fi

(root -b -l  <<EOF
.L filter.C
filter($pid,"${sub_dir}/${root_list[$i]}", "${sub_dir}/${four_vec}_${i}.root")
EOF
)&
done
wait

hadd ${sub_dir}/"$four_vec".root ${root_four_vec[@]}

#Be careful with the remove (rm) command!
rm ${sub_dir}/${four_vec}_*.root
}



if [[ $pythia = on ]]; 	  then pythia;    fi
if [[ $lhe2root = on ]];  then lhe2root;  fi
if [[ $p4_filter = on ]]; then p4_filter; fi
