# OpenMP thread-affinity
export OMP_PROC_BIND=spread
export OMP_PLACES=threads

# general variable declarations
tstep=0.002
inpxlength=30
inpylength=30
inpzlength=30
inpnumAtoms=9990
temp=1.5

# thermodynamic force variable declarations
densbinwidth=0.2
neighbors=0
damping=0.0
forcemod=1.0

# region management 
centralmin=0.0
centralmax=0.0
cappingmin=0.0
cappingmax=12.5
thermostatmin=0.0
thermostatmax=15.0
thermoforcemin=10.0
thermoforcemax=14.5

# manage input 
inpfile=equilibrateLangevin_final.h5md

# manage output
outfile=thermoForce_naive_n"$inpnumAtoms"_T"${temp//./}"_$(date "+%Y_%m_%d")
outparent=rho0370_T15_capAll_thermAll_2026_04_01
outdir=/srv/public/julianhille/project/mrmd_simulations/$outparent/$outfile

# non-critical metadata
echo "hostname: $(hostname)"
echo "date: $(date)"
echo "working directory: $(pwd)"
echo "tstep: $tstep"
echo "temp: $temp"
echo "densbinwidth: $densbinwidth"
echo "neighbors: $neighbors"
echo "damping: $damping"
echo "forcemod: $forcemod"
echo "centralmin: $centralmin"
echo "centralmax: $centralmax"
echo "cappingmin: $cappingmin"
echo "cappingmax: $cappingmax"
echo "thermostatmin: $thermostatmin"
echo "thermostatmax: $thermostatmax"
echo "thermoforcemin: $thermoforcemin"
echo "thermoforcemax: $thermoforcemax"
echo "inpfile: $inpfile"
echo "forceguess: $forceguess"
echo "outfile: $outfile"
echo "outdir: $outdir"

mkdir $outdir
cp $inpfile $outdir
cp $0 "$outdir"/run_script.sh

# thermodynamic force
./04_LennardJones_GentleChemostat --inpfile $inpfile --outint 5000000 --outfile "$outdir"/"$outfile" --nsteps 200000001 --sampling 500 --update 5000000 --centralmin $centralmin --centralmax $centralmax --cappingmin $cappingmin --cappingmax $cappingmax --thermostatmin $thermostatmin --thermostatmax $thermostatmax --thermoforcemin $thermoforcemin --thermoforcemax $thermoforcemax --tstep $tstep --temp $temp --densbinwidth $densbinwidth --neighbors $neighbors --damping $damping --forcemod $forcemod > "$outdir"/"$outfile".out