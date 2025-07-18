pwd=$PWD
Simu_DIR=/home/llr/ilc/shi/code/SiWECAL-Prototype/generation/run_scripts/Simu2025-06
#source /cvmfs/sw.hsf.org/key4hep/setup.sh -r 2025-05-29
. /cvmfs/sw.hsf.org/key4hep/setup.sh
export PATH=/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/python-venv/1.0-2b64yh/bin:$PATH

echo "flatpak version"
FAKE_LDAP_PATH="${Simu_DIR}/lib/fake_ldap.so"
if [ ! -f "$FAKE_LDAP_PATH" ]; then
    echo "Error: fake_ldap.so not found at $FAKE_LDAP_PATH"
    exit 1
fi
export LD_LIBRARY_PATH=${FAKE_LDAP_PATH}:$LD_LIBRARY_PATH
ldd $(which flatpak) | grep libldap

#. source /home/llr/ilc/shi/key4hep/env.sh
echo "root version"
which root
echo "python version"
which python
python --version

additional_paths=(
  "/home/llr/ilc/shi/code/SiWECAL-Prototype/generation/run_scripts/Simu2025-06/build/lib"
  "/home/llr/ilc/shi/code/SiWECAL-Prototype/generation/run_scripts/Simu2025-06/build/lib64"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/libtool/2.4.7-cgp2sf/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/root/6.32.04-vms5ij/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/vdt/0.4.4-feme7a/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/root/6.32.04-vms5ij/lib/root"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/gsl/2.7.1-v2sxnu/lib/root"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/opendatadetector/v3.0.0-sdtl34/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/opendatadetector/v3.0.0-sdtl34/lib64"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4rectracker/0.5.0-w3rrqg/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4reco/0.2.1-aaz2ba/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4mljettagger/0.1.0-zfjkox/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4gaudipandora/0.1.0-vp7mlu/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4clue/01-00-07-b2o4iz/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4marlinwrapper/00-12-hxzsgw/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4simdelphes/00-07-05-5yyzax/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4reccalorimeter/0.1.0pre17-o3kzzx/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4simgeant4/0.1.0pre16-sp5khc/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4geo/00-22-ubhvqv/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/lcio/2.22.6-zicklo/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/fccanalyses/0.11.0-a25snu/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/fccdetectors/0.1pre10-hs76r7/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/edm4hep/0.99.2-sgaasm/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/podio/1.3-njsv57/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4gen/0.1pre14-iyq5vx/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/k4fwcore/1.3-lix236/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2025-05-29/x86_64-almalinux9-gcc14.2.0-opt/dd4hep/1.32-x2e5mf/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/sio/0.2-lp6t6f/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/r/4.4.1-32pivc/rlib/R/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/hdf5/1.14.3-pgt77p/lib"
  "/cvmfs/sw.hsf.org/contrib/x86_64-almalinux9-gcc11.4.1-opt/gcc/14.2.0-yuyjov/lib64"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/gosam-contrib/2.0-xeevlf/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/intel-oneapi-mkl/2024.2.1-coyg5v/mkl/2024.2/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/zlib-ng/2.2.1-g7eiqc/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/xz/5.4.6-qpg4p7/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/xxhash/0.8.2-qbue5z/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/lz4/1.9.4-cozw3x/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/zstd/1.5.6-h3ocxj/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/intel-tbb/2021.12.0-utctht/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/ncurses/6.5-z23h5j/lib"
  "/cvmfs/sw.hsf.org/key4hep/releases/2024-10-03/x86_64-almalinux9-gcc14.2.0-opt/openssl/3.3.1-ju4slh/lib64"
)


for p in "${additional_paths[@]}"; do
    if [[ ":$LD_LIBRARY_PATH:" != *":$p:"* ]]; then
        export LD_LIBRARY_PATH="$LD_LIBRARY_PATH:$p"
        echo "Added $p to LD_LIBRARY_PATH"
    else
        #echo "Path $p already in LD_LIBRARY_PATH"
        :
    fi
done

cd $Simu_DIR/build
echo "k4 local repo"
k4_local_repo .
cd $pwd
