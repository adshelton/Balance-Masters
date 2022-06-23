# Balance-Masters
Master set of codes used in balance perturbation research (Currently set for Per Motor Rep Project 2022)
Requires anc and trc files

## CURRENT LIST OF CODES & FUNCTION
### Main codes:
UNCABL_convert_anc -> Convert anc files from Cortex using .cal file

UNC_Balance_master_run1 -> Creates base matalb structure for subject (current version is for per motor rep session 1)

UNC_Balance_master_run2 -> Balance metrics (work in progress)

UNC_balance_master_run3_lle -> short and long term divergence exponents for local dynamic instability

MoS_AllSubjects_v3 -> Stand alone code that does not require run1 to provide margin of stability at for trials at points of importance

Slip_WBAM_prog -> Whole body angular momentum rnages for treadmill slips. Shows progression and difference from 1st to 5th slip

### Auxillary Codes:
slip_occ -> produces matlab structure for HS that triggers slip perturbation

LFP_occur_left & LFP_occur_right -> produces matlab structure for HS that starts gait cycles relevant to waist pull

### Necessary Functions:
lyarosenstein -> function needed for UNC_balance_master_run3_lle

load_zero

load_trc

load_anc

load_fpcal

convertFPdata
