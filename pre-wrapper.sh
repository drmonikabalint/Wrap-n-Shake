#!/bin/bash
			#
			#	
			# PRE-WRAPPER 1.1
			#	
			# Preparing process, Pre-Wrapper Stage
			#			
			# If use wrapper.sh, please cite: 
			#
			#                    Balint M. et. al, (2017) Journal of chemoinformatics 
			#
			# USAGE: ./pre-wrapper.sh -t TargetName(.pdb)
			#		          -l LigandName(.pdb)
			#			  -p parameter path
			#		      
			# 		
			# M. Bálint, 2017

echo " 			PRE-WRAPPER Version 1.1,  19.12.2017"

#		
# Assuming that moving forward, you have a minimized ligand and target, with all the hydrogen added, corresponding to the 
# total charge of the molecule, this script will generate the:
#				1. pdbqt files (ligand and target) and 
#				2. parameter files (grid and docking) 
#


while [[ $# > 1  ]]
do
	key="$1"
	case $key in
    		-t|--target)
    		target="$2"
		
		###
		### This if condition is essential, to declare a new variable, called "target_name", with the pdbqt extension cut off.
		### "target_name" variable will be used in the cycle function, when renaming the target files.  	
		###
		
		if [[ "$target" ]]; then
		target_name=${target%\.pdb}
		fi

    		shift                       	# past argument
    		;;

		-l|--ligand)
		ligand="$2"

		if [[ "$ligand" ]]; then
		ligand_name=${ligand%\.pdb}
		fi

		shift 				# past argument
		;;


 		-p|--SCRIPTPATH)
    		SCRIPTPATH="$2"
    		shift                           # past argument
    		;; 

    		*)
	
	esac 
	shift                            	# past argument or value

done


echo Target name   = "${target}"
echo Ligand name   = "${ligand}"
echo Scripts and parameterpath = $PARAMETERPATH


if [[ -z $target ]] || [[ -z $ligand ]] || [[ -z $SCRIPTPATH ]]
then

	echo "***ERROR***"
	echo "***Incomplete arguments provided!" 
	echo "***Correct usage of the script: $0 -t Target Name -l Ligand Name -p script and parameter path"

	exit 1
fi

	# run the python scripts with customized parameters.

	# prepare receptor
	$SCRIPTPATH/pythonsh $SCRIPTPATH/prepare_receptor4.py -r $target -o $target_name.pdbqt   -v -U nphs
	# prepare ligand
	$SCRIPTPATH/pythonsh $SCRIPTPATH/prepare_ligand4.py -l $ligand  -o $ligand_name.pdbqt -v -d  $SCRIPTPATH/ligand_dict.py -F 
	# prepare docking parameter file
	$SCRIPTPATH/pythonsh $SCRIPTPATH/prepare_dpf42.py 	-l $ligand_name.pdbqt -r $target_name.pdbqt -p ga_num_evals=20000000 -p ga_pop_size=250 -p ga_num_generations=2000000000 -p ga_run=100 -p rmstol=2
	# prepare grid parameter file
	$SCRIPTPATH/pythonsh $SCRIPTPATH/prepare_gpf4.py 	-l $ligand_name.pdbqt -r $target_name.pdbqt -p spacing=0.375 -p npts='200,200,200' -p ligand_types='A,C,F,H,HD,HS,N,NA,NS,OA,OS,S,SA,P,MG,MN,Zn,CA,BR,Cl,YY,LL' -v
	#
	# Some changes needs to be done on the prepared DPF and GPF files.
	#
	mv $ligand_name\_$target_name.dpf $target_name.dpf	
	############
	# Add the new parameter file name, in the first row. This is how the LL and YY atom parameters will be read by AutoDock4, instead of the default parameter file.
	sed -i -e "1 c\parameter_file $SCRIPTPATH/AD4_parameters.dat     #new atom types defined"  $target_name.dpf
	# Add more atom types. The same atom types, that were added in the gpf file.
	sed -i -e "s/ligand_types.*$/ligand_types A C F H HD HS N NA NS OA OS S SA P MG MN Zn CA BR Cl YY LL     # receptor atom types /g"  $target_name.dpf
	sed -i -e "/^map $target_name/"d $target_name.dpf
	sed -i -e "7i map $target_name.A.map                 # atom-specific affinity map"   $target_name.dpf
	sed -i -e "8i map $target_name.C.map                 # atom-specific affinity map"   $target_name.dpf
	sed -i -e "9i map $target_name.F.map                 # atom-specific affinity map"   $target_name.dpf
	sed -i -e "10i map $target_name.H.map                # atom-specific affinity map"   $target_name.dpf
	sed -i -e "11i map $target_name.HD.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "12i map $target_name.HS.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "13i map $target_name.N.map                # atom-specific affinity map"   $target_name.dpf
	sed -i -e "14i map $target_name.NA.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "15i map $target_name.NS.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "16i map $target_name.OA.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "17i map $target_name.OS.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "18i map $target_name.S.map                # atom-specific affinity map"   $target_name.dpf
	sed -i -e "19i map $target_name.SA.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "20i map $target_name.P.map                # atom-specific affinity map"   $target_name.dpf
	sed -i -e "21i map $target_name.MG.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "22i map $target_name.MN.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "23i map $target_name.Zn.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "24i map $target_name.CA.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "25i map $target_name.BR.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "26i map $target_name.Cl.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "27i map $target_name.YY.map               # atom-specific affinity map"   $target_name.dpf
	sed -i -e "28i map $target_name.LL.map               # atom-specific affinity map"   $target_name.dpf
	#############
	# Add the new parameter file name, in the first row. This is how the LL and YY atom parameters will be read by AutoDock4, instead of the default parameter file.	
	sed -i -e "1i parameter_file $SCRIPTPATH/AD4_parameters.dat     #new atom types defined" *.gpf

	# Add more atom types, so that this gpf file can be used in the more advanced wrapper cycles, when additional atoms types will also be present on the protein 
	sed -i -e 's/receptor_types.*$/receptor_types A C F H HD HS N NA NS OA OS S SA P MG MN Zn CA BR Cl YY LL      # receptor atom types /g' *.gpf

