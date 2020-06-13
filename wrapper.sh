#!/bin/bash

			#	
			# WRAPPER 1.1
			# If use wrapper.sh, please cite: 
			#
			#                    Balint M. et. al, (2017) Journal of chemoinformatics  
			#
			# USAGE: ./wrapper.sh -n CycleNumber
			#		      -t TargetName(.pdbqt)
			#		      -l LigandName(.pdbqt)
			#		      -b Binary path (ex: /home/user/bin/)
			#		      -p Parameter path (ex: /home/user/bin/)		
			# 		      -r ligand template(.pdbqt)
			# Bálint M, 2017


echo " 			WRAPPER Version 1.1,  02.01.2018"

####################################################################################################################
#
# Declare the variable from the user provided arguments
# 
####################################################################################################################

end_cycle=10000

while [[ $# > 1  ]]
do
	key="$1"
	case $key in
    		-n|--end_cycle)
    		end_cycle="$2"
    		shift                           # optional argument
    		;;

    		-t|--target)
    		target="$2"
		
		###
		### This if condition is essential, to declare a new variable, called "target_name", with the pdbqt extention cut off.
		### "target_name" variable will be used in the cycle function, when renaming the target files.  	
		###
		
		if [[ "$target" ]]; then
		target_name=${target%\.pdbqt}
		fi

    		shift                       	# past argument
    		;;

		-l|--ligand)
		ligand="$2"

		###
		### This if condition is essential, to declare a new variable, called "ligand_name", with the pdbqt extention cut off.
		### "ligand_name" variable will be used in the cycle function, when renaming the target files.  	
		###

		if [[ "$ligand" ]]; then
		ligand_name=${ligand%\.pdbqt}
		fi
		shift                       	# past argument
    		;;


    		-b|--BINARYPATH)
    		BINARYPATH="$2"
    		shift                           # past argument
    		;; 

   		-p|--PARAMETERPATH)
    		PARAMETERPATH="$2"
    		shift                           # past argument
    		;; 

   		-r|--template)
    		template="$2"
    		shift                           # past argument
    		;; 


    		*)
	
	esac 
	shift                            	# past argument or value

done

echo Cycle number  = "${end_cycle}"
echo Target name   = "${target}"
echo Ligand name   = "${ligand}"
echo Parameter path = "$PARAMETERPATH"
echo Binary path    = "$BINARYPATH"
echo Template = "${template}"


####################################################################################################################
# Function named "cycle_docking" is created, in order to perform the same steps in multiple docking cycles 
# Cycle function has 3 arguments. 
#                      -1st: number of the current cycle
#		       -2nd: the name of the current target 
#		       -3rd: number of the next cycle
####################################################################################################################

function cycle_docking
{

			# STEP 1
			# in each cycle docking process is performed in the first step
			#

			echo "******"	
			echo "****** Wrapper process continues with cycle number: $1!"		
			echo "****** EXIT Condition was not yet reached "
			echo "****** Started AutoGrid"	

			$BINARYPATH/autogrid4  -p $PARAMETERPATH/$target_name.gpf -l $target_name.glg
			echo "****** Started AutoDock"	
			$BINARYPATH/autodock4  -p $PARAMETERPATH/$target_name.dpf -l $target_name.dlg  
			mv $target_name.glg  $2\_wrp.glg
			mv $target_name.dlg  $2.dlg
			mv $target_name.pdbqt  $2\_wrp.pdbqt	
			rm *map*							
			
			# STEP 2
			# Cluster generation, ranking, and atom type assignation (YY and LL), using wrp. 

			$BINARYPATH/wrp -f $2.dlg  -difc 3.5 -drnk 2 -dsgn 3.5  -v silent -t $2\_wrp.pdbqt -c $1	
			mv $2.dlg $2\_wrp.dlg
}

####################################################################################################################
# Function named "cycle_eval" is created, in order 
# to calculate the accessible surface area 
# and
# to evaluate the outputs of the docking cycle
####################################################################################################################


function cycle_eval
{
			# STEP 3
			echo "****** Surface calculations"
gmx editconf -f O_$2\_input_t.pdb -o O_$2\_input_t.gro -d 1 -bt cubic > out1.txt 2>&1
gmx editconf -f O_$2\_input_t.gro -o O_$2\_input_t_sas.pdb > out2.txt 2>&1
gmx sasa -f O_$2\_input_t_sas.pdb -s O_$2\_input_t.pdb  -oa total_atomarea.xvg -surface -probe 0.35 > out3.txt 2>&1 <<EOF
0
EOF

			rm out1.txt
			rm out2.txt
			rm out3.txt

			awk '{print $4}' O_$2\_input_t.pdb > res_name.txt
			sed  -i -e '/^$/d' res_name.txt
			sed  -i -e '/#/d' total_atomarea.xvg
			sed  -i -e '/@/d' total_atomarea.xvg
			rm O_$2\_input_t.gro
			rm O_$2\_input_t_sas.pdb

		        paste res_name.txt total_atomarea.xvg >  total_atomarea_lig.xvg

			# from the obtained xvg file, we remove the lines that are containing the surface of the ligands, keeping the surface calculated strictly for the receptor residues

			sed '/LIG/,$d' total_atomarea_lig.xvg > O_$2\_ASA.log

			# calculate the sum of the ASA 
			awk '{print  sum1+=($3)}'  O_$2\_ASA.log > tmp
			cp tmp O_$2\_ASA.log
			sed '$!d' O_$2\_ASA.log > O_$2\_ASA_free.log		

			# STEP 4
			# the concatenated receptor is then moved to the next cycle serving as the receptor inpus for docking.

			cp O_${target_name}_${i}_wrp.pdbqt $target_name.pdbqt	
			mv $target_name.pdbqt ../$3
			
			rm *ASA.log 
			rm tmp
			rm total_atomarea.xvg
			rm total_atomarea_lig.xvg
			rm res_name.txt
			rm area.xvg



			#Step 5 compare the ASA. 
			#Perform the lines in the if condition, only if the folder is not equal to 1. 
			#This is a comparison, and calculation of how much (expressed in percentage) accessible surface area remained on the target, compared to the reference (target in round 1, without any ligand)
			
			if [[ "${PWD##*/}" != "1" ]];then

				cp ../1/O_${target_name}_1_ASA_free.log  ASA_reference.log 	
				paste  ASA_reference.log O_$2\_ASA_free.log | awk '{print ($2*100)/$1}' > O_$2\_surface_percentage.log
				awk '{if ($1 <= 1) print $1}'  O_$2\_surface_percentage.log  > O_$2\_surface_exit.log
				
				echo "Available free protein surface after cycle number $i:" 
				cat O_$2\_surface_percentage.log 
		
				mv  O_$2\_surface_percentage.log ../stats
				rm ASA_reference.log
				

			fi

			#Step 6 check the energy
			echo "****** Energy evaluations"
			sed -n 2p O_${target_name}_${i}_wrp.sta | awk '{if ($3 >= 0) print $3}' > O_$2\_energy_exit.log

			sed -n 2p O_${target_name}_${i}_wrp.sta > tmp
			echo "Lowest energy after cycle number $i:" 
			awk '{print $3}' tmp

			awk '{print $1, "\t", $3, "\t" }' tmp > O_$2\_lowest_energy.log		
			rm tmp
			mv O_$2\_lowest_energy.log ../stats
			

			find . -size 0 -delete
}


####################################################################################################################
#
# Variables declared
# The following if condition, is a safety check. Verifies if user provided all mandatory arguments.
#
####################################################################################################################

if [[ -z $target ]] || [[ -z $ligand ]] || [[ -z $BINARYPATH ]] || [[ -z $PARAMETERPATH ]] || [[ -z $template ]]
then

	echo "***ERROR***"
	echo "***Incomplete arguments provided!" 
	echo "***Correct usage of the script: $0 -n Cycle Number (optional) -t Target Name -l Ligand Name -b binary path -p parameter path -r template "

	exit 1
fi


####################################################################################################################
# Preparatory checks and variable definitions were done. 
# The next part of the script does the actual work.
# 
# At start, we create directories, where the docking cycles will take plase, in order to have a continuous process
#
# All the necessary inputs are copied  to the current working folders
#
# NOTE 1 please note and remember, that TargetName.pdbqt is only copied to directory nr. 1 since this is the only cycle
# 	 where the protein receptor will be used in this form. For the following cycles, target will have
# 	 some of the atom types renamed ( see atom YY ) and will have the ligands from previous cycles concatenated.
#
# NOTE 2 please note that *ref.pdb (reference file) is only copied into the folders when RMSD calculation is needed 
#	 to copy the reference molecule, please remove the hashtag "#" from the below for loop.  
#	 to calculate the RMSD please add the -r argument and the reference file name to wrp command line in the function cycle.  
#
####################################################################################################################





####################################################################################################################
##
## Functions "cycle_docking" and "cycle_eval" are called from here. 
##
## Exit criterions are checked from an if condition inside the while loop.
##
####################################################################################################################

	# loop will terminate, if one of the two conditions are met.
	# condition 1 is the surface exit
	# condition 2 is the energy exit 
i=0
mkdir 1
mkdir stats
while [[ "$i" != $end_cycle ]] &&  [[ ! -f *_surface_exit.log ]] &&  [[ ! -f *_energy_exit.log ]]; 
do
	(( i++ ))	
	
	next=`expr "$i" + 1`;

	mkdir "${next}"
#       cp *ref.pdb ${i}/
	cp $target_name.pdbqt 1/
	cp $ligand_name.pdbqt ${i}/

	cd ${i}

  	outfile="O_${target_name}_${i}_wrp.pdbqt"
        if [ -e $outfile ]; then
                        echo "******"  
                        echo "****** Cycle $i is already finished: SKIPPING !"               
                        echo "****** Reason: found file $outfile  "
                        echo "******"  
			cycle_eval "${i}"  "${target_name}_${i}" "${next}"
			
        else
                        cycle_docking "${i}"  "${target_name}_${i}" 
			cycle_eval "${i}"  "${target_name}_${i}" "${next}"
			
        fi

	if  [ -e *_surface_exit.log ]; then              

		echo "******"	
		echo "****** Wrapper process is soon to be finished!"		
		echo "****** SURFACE EXIT Condition was reached at the end of Cycle number: $i"
		echo "******"	
		
		# wrp will calculate the distance between the cluster representatives and the target, 
		# wrp will eliminate ligands that are above -dmax cut-off 

		cp ../$template  ../${i}

		$BINARYPATH/wrp -f O_${target_name}_${i}\_wrp.pdbqt  -dmax 3.5 -m trimming -p $template

		mv O_O_${target_name}_${i}\_wrp_trm.pdb O_${target_name}_${i}\_wrp_trm.pdb
		mv O_O_${target_name}_${i}\_wrp.log O_${target_name}_${i}\_wrp.log
		cd ..
	
		break
	fi

	if [  -e *_energy_exit.log ] ; then              

		echo "******"	
		echo "****** Wrapper process is soon to be finished!"		
		echo "****** ENERGY EXIT Condition was reached at the end of Cycle number: $i"
		echo "******"			

		# wrp will calculate the distance between the cluster representatives and the target, 
		# wrp will eliminate ligands that are above -dmax cut-off 

		cp ../$template  ../${i}

		$BINARYPATH/wrp -f ${target_name}_${i}\_wrp.pdbqt  -dmax 3.5 -m trimming -p $template
		
		cd ..
		break
	fi

	
	cd ..

done

cd stats 
cat *_lowest_energy.log > lowest_energy_per_cycle.log
cp lowest_energy_per_cycle.log ../stats
cd ../stats
echo "100.0000" > O_${target_name}_1_surface_percentage.log
cat *_surface_percentage.log > surface_percentage_per_cycle.log
echo "Cycle  Energy		ASA"
echo "Nr.    (kcal/mol)	(%)"

paste lowest_energy_per_cycle.log  surface_percentage_per_cycle.log 

cd ..
