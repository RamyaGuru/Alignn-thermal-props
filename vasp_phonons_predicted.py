#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  9 13:18:40 2022

@author: rlg3


Set up predicted compound phonon calculations
"""

from jarvis.tasks.vasp.vasp import (
    JobFactory,
    VaspJob,
    GenericIncars,
    write_jobfact,
)
from jarvis.io.vasp.inputs import Potcar, Incar, Poscar
from jarvis.db.jsonutils import dumpjson
from jarvis.db.figshare import data
from jarvis.core.atoms import Atoms
from jarvis.tasks.queue_jobs import Queue
import os

##############################################################################
# Make sure you have latest version of jarvis-tools, pip install -U jarvis-tools
# Specify your vasp_cmdcluster type etc. info here
# VASP_PSP_DIR should be defined in the PATH as pseudopotential directory

vasp_cmd = "/usr/bin/mpirun vasp_std"

# Change to your path of .bindat file
copy_files = ["/users/rlg3/bin/vdw_kernel.bindat"]

# For slurm
submit_cmd = ["sbatch", "submit_job"]
##############################################################################
def get_atoms(jid="", mpid="", oqmd_id="", aflow_id=""):
    """
    Provide only one of these IDs.
    Examples:
    jid='JVASP-1002' or mpid='mp-149' or
    oqmd_id='10215', or aflow_id='48708a4622918820'
    """
    if mpid != "":
        mp = data("mp_3d")
        for i in mp:
            if i["id"] == mpid:
                atoms = Atoms.from_dict(i["atoms"])
                del mp
                return atoms
    if jid != "":
        jv = data("dft_3d")
        for i in jv:
            if i["jid"] == jid:
                atoms = Atoms.from_dict(i["atoms"])
                del jv
                return atoms
    if oqmd_id != "":
        oq = data("oqmd_3d")
        for i in oq:
            if i["id"] == oqmd_id:
                atoms = Atoms.from_dict(i["atoms"])
                del oq
                return atoms
    if aflow_id != "":
        af1 = data("aflow1")
        for i in af1:
            if i["id"] == aflow_id:
                atoms = Atoms.from_dict(i["atoms"])
                del af1
                return atoms
        af2 = data("aflow2")
        for i in af2:
            if i["id"] == aflow_id:
                atoms = Atoms.from_dict(i["atoms"])
                del af2
                return atoms


# If a user wants to run on its on ,aterials
# atoms = Poscar.from_file('YourPOSCAR')
# and send it to atoms in the script below
# Select/desect you want to run
# More detais in
# https://github.com/usnistgov/jarvis/blob/master/jarvis/tasks/vasp/vasp.py#L81
steps = [
    "ENCUT",
    "KPLEN",
#    "RELAX",
#    "BANDSTRUCT",
#    "LOPTICS",
#    "MBJOPTICS",
    "ELASTIC",
]
incs = GenericIncars().optb88vdw().incar.to_dict()


# List of materials divided into chunks of 50


ids = ["JVASP-112289", "JVASP-25379"]


for id in ids:
    atoms = get_atoms(jid=id)
    mat = Poscar(atoms)
    mat.comment = "bulk@" + str(id)
    cwd_home = os.getcwd()
    dir_name = id + "_" + str("PBEBO")
    if not os.path.exists(dir_name):
        os.makedirs(dir_name)
    os.chdir(dir_name)
    job = JobFactory(
        vasp_cmd=vasp_cmd,
        poscar=mat,
        steps=steps,
        copy_files=copy_files,
        use_incar_dict=incs,
    )

    dumpjson(data=job.to_dict(), filename="job_fact.json")
    write_jobfact(
        pyname="job_fact.py",
        job_json="job_fact.json",
        input_arg="v.step_flow()",
    )

    # # Example job commands, need to change based on your cluster
    # job_line = (
    #     "source ~/rk2/rlg3/mambaforge/envs/my_jarvis/bin/activate my_jarvis \n"
    #     + "python job_fact.py"
    # )
    # name = id
    # directory = os.getcwd()
    # # Queue.pbs(
    # #     job_line=job_line,
    # #     jobname=name,
    # #     directory=directory,
    # #     submit_cmd=submit_cmd,
    # # )
    # # os.chdir(cwd_home)


    # # For Slurm clusters
    # Queue.slurm(
    #     walltime = '10:00:00',
    #     queue='main',
    #     job_line=job_line,
    #     jobname=name,
    #     jobout = name + '.out',
    #     joberr = name + '.err',
    #     directory=directory,
    #     submit_cmd=submit_cmd,
    # )
    # os.chdir(cwd_home)

