# Panda Analysis

Throughout this readme, I will refer to several user-defined environment variables. 
These are typically defined in `T3/setup.sh`, but the user can define things elswhere.

## Installation

```bash
cmsrel CMSSW_9_3_0
cd CMSSW_9_3_0_patch1/src
cmsenv
git clone https://github.com/PandaPhysics/PandaTree
git clone https://github.com/sidnarayanan/PandaCore
git lfs clone https://github.com/sidnarayanan/PandaAnalysis
PandaCore/bin/genDict -f -j8                                          # I typically add PandaCore/bin to my $PATH
scram b -j8
```

Explanation: `PandaTree` is the data format, `PandaCore` is some core utilities for dealing with ROOT and python, and `PandaAnalysis` implements the analyses.
Most of your interaction should only be with `PandaAnalysis`, unless you find a bug, which you will.

## Producing a flat tree for analysis

There is a core analysis tool (`PandaAnalyzer`) that outputs a data format (`GeneralTree`).
`PandaAnalyzer` is configured using an `Analysis` class, which can turn on and off analysis-specific calculations and corresponding branches in the output tree.
The functions are implemented in `PandaAnalysis/Flat/src/Modules*cc`.

For what follows, I'll assume you're inside `$CMSSW_BASE/src/PandaAnalysis`.

To see what is produced in the output tree, you can look in `Flat/config/GeneralTree.cfg`. 
If you want to add a variable, just put it in the config and then run:
```bash
./config/generateTreeClass.py --config config/GeneralTree.cfg
```
Any custom code in the class definition will be preserved, and the new variables will be added.
By default, the new variables are booked always, but if you want to make it conditional (e.g. only if the VBF flag is set), you can modify by hand `GeneralTree::Write`.

### Defining your analysis

The analysis framework is heavily integrated with the MIT T3/T2. 
In princple it can also run on SubMIT, but it gets tricky (submission is on one node, stageout is on another).
Until the T3 drives are visible from submit.mit.edu, or I figure out remote condor submission (still debugging this...), let's assume you only run on the T3/T2.
For most vanilla analyses, this is fine, as a typical analysis takes ~2 hours, and the failure rate on the T3 is essentially zero (much higher on SubMIT).
To run only on the T3, set `export SUBMIT_CONFIG=T3`; to include the T2, do `export SUBMIT_CONFIG=T2`.

For what follows, I assume you're in `$CMSSW_BASE/src/PandaAnalysis/T3`.

First, open up `setup.sh`. 
This defines all the environment variables we'll need.
Make sure anything that is a path (`PANDA_FLATDIR`,`SUBMIT_LOGDIR`,etc) is writable by you.
The `SUBMIT*` environment variables are used for running the jobs, and the `PANDA*` variables are for running things on the outputs of the jobs.
`SUBMIT_TMPL` points to a script inside `inputs/`, that will define your analysis.
For a straightforward example, do `export SUBMIT_TMPL=skim_gghbb_tmpl.py`.

Now, let's open up `inputs/skim_wlnhbb_tmpl.py` as a concrete example and look at it.
The main function you have to worry about is `fn`.
The rest is all fluff.
For convenience, here's the (annotated) content of `fn`:
```python
def fn(input_name, isData, full_path):

    PInfo(sname+'.fn','Starting to process '+input_name)
    # now we instantiate and configure the analyzer
    skimmer = root.PandaAnalyzer()
    analysis = wlnhbb(True)                                          # this is imported from PandaAnalysis.Flat.Analysis, where the defaults are set
    analysis.processType = utils.classify_sample(full_path, isData)  # set the type of the process
    if analysis.processType == root.kTT or analysis.processType == root.kSignal:
        analysis.reclusterGen = True                                 # additional customization can be done on the fly to the analysis object
    skimmer.SetAnalysis(analysis)
    skimmer.isData=isData
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kVHBB)             # set the preselection
    skimmer.SetPreselectionBit(root.PandaAnalyzer.kPassTrig)         # only save data events that trip a trigger

    return utils.run_PandaAnalyzer(skimmer, isData, input_name)      # run the analysis 
```

### Testing your analyzer

Inside `Flat/test`, there is a testing script that runs as:
```bash
./test.py /path/to/input/panda.root [DEBUG_LEVEL]
```
You have to open up the script and modify the number of events you want to run, the type of file it is (data, W+jets MC, etc), and what flags are on.

## Running on the grid

### Cataloging inputs

```bash
./catalogT2Prod.py --outfile /path/to/config.cfg [ --include datasets to include ] [ --exclude datasets to skip ] [ --force ] [--smartcache]
```

I recommend you put the config in a web-accessible place for use in later steps. For example:
```bash
./catalogT2Prod.py --force --outfile ~/public_html/histcatalog/$(date +%Y%m%d).cfg --include TT --exclude TTbarDM --smartcache
```

The above command will do the following things:

- It will only check datasets that contain `TT` and do not contain `TTbarDM` in the dataset's nickname

- If a dataset is not found in `PandaCore.Tools.process`, it will guess the nickname of the dataset, give it a xsec of 1, and write it to the catalog (`--force`)

- If the file does not exist locally on the T3, a smartcache request will be made

- The output will be a timestamped web-facing config file


### Building the work environment

First, make sure the following environment variables are defined (some examples shown):
```bash
export PANDA_CFG="http://snarayan.web.cern.ch/snarayan/eoscatalog/20170127.cfg"  # location of config file from previous section
export SUBMIT_TMPL="skim_merge_tmpl.py"  # name of template script in T3/inputs
export SUBMIT_NAME="v_8024_2_0"  # name for this job
export SUBMIT_WORKDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/work/"  # staging area for submission
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/logs/"  # log directory
export SUBMIT_LOGDIR="/data/t3serv014/snarayan/condor/"${SUBMIT_NAME}"/locks/"  # lock directory
export SUBMIT_OUTDIR="/mnt/hadoop/scratch/snarayan/panda/"${SUBMIT_NAME}"/batch/"  # location of unmerged files
export PANDA_FLATDIR="${HOME}/home000/store/panda/v_8024_2_0/"   # merged output
export SUBMIT_CONFIG=T2  # allow running on T3 or T2. if $SUBMIT_CONFIG==T3, then only run on T3
```

`T3/inputs/$SUBMIT_TMPL` should be the skimming configuration you wish to run your files through. 

### Testing the jobs

If you want, you can test-run a job locally:
```bash
./task.py --build_only --nfiles 1   # just build the job, with one file per job
cd $SUBMIT_WORKDIR                  # go to the working directory
python skim.py 0 0                  # run the first job in the configuration
cd -                                # back to bin
./task.py --clean                   # clean up after yourself
```

### Submitting and re-submitting

To submit jobs, simply do
```bash 
./task.py --submit [--nfiles NFILES] [--clean]
```
where NFILES is the number of files in each job. 
The default is 25 files if that flag is not passed.
The clean argument will make sure to wipe out all staging directories to ensure a clean release (it's optional because sometimes you don't want to do this).

To check the status of your jobs, simply do:
```bash
./task.py --check [--silent] [--force] [--nfiles NFILES] [--monitor NSECONDS]
```
Note that the above command overwrites `$SUBMIT_WORKDIR/local.cfg`, with the intention of preparing it for resubmission.
The file will be recreated as a configuration to rerun files that are not present in the output and not running.
- `--silent` will skip the per-sample breakdown.
- `--nfiles` will repackage `local.cfg` into a different number of files per job.
- `--force` will re-catalog files that are incomplete, not just missing.
- `--monitor` will capture the terminal screen and refresh the status if either (a) a job has completed succesfully or (b) `NSECONDS` has elapsed since the last refresh.

To resubmit missing files, simply do
```bash
./task.py --silent
```
In the case that you are using the `--force` option, make sure you have no running jobs before resubmitting, or you may end up with duplicated outputs.
Forcing resubmission is generally discouraged and is largely included for historical reasons.
The job framework is robust enough at this point that forcing should not be necessary.

If you are having lots of failures, it may be interesting to analyze the logs located in `$SUBMIT_LOGDIR`. 
You can do this manually, or:
```bash
./analyzeLogs.py [--dump]
```
This will print to screen a basic analysis of the errors observed, which errors are correlated, and how frequently they occur.
If the `--dump` flag is passed, then a directory `log_dumps/` is created, containing detailed information on each failure class (where it failed and on for what inputs).


## Merging

Make sure `$PANDA_FLATDIR` exists. Then, go into `T3/merging` and do:
```bash
./merge.py [--cfg CONFIG] [--silent] TTbar_Powheg
```
to merge the Powheg TT sample, for example. 
If provided, `CONFIG` is the module that is imported from `configs/`. 
The default is `common`, but there are others, like `leptonic`.
To merge en-masse (e.g. many many signal outputs), you can do something like:
```bash
submit --exec merge.py --arglist list_of_signals.txt
```
This assumes that you have `PandaCore/bin` in your `$PATH`
The `submit` command will print a cache directory it created to keep track of the jobs.
You can check the status of the jobs by doing
```bash
check --cache <cache_directory> [--resubmit_failed]
```
The last flag is optional and will resubmit anything that failed (exited without code 0).
NB: the `submit` and `check` executables are very generic and don't know anything about the code they are running.
For them, "success" simply means "exited with code 0".
So it is important to check that the output looks sane to you.
