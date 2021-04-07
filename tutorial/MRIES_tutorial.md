# Welcome to the MRIES tutorial
This tutorial is an introduction to MRIES processing and visualization.

The MRIES Tutorial is split into five parts, installing MRIES, data inputting, processing pipline, visualization and result outputting.


## 1. Installing MRIES:
1. Start Matlab.

2. Add the *MRIES* folder to path of Matlab. 

<img  src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/path.png"/>

3. Type "MRIES" at the Matlab command prompt ">>" and press enter. If MRIES main GUI will be shown, MRIES will be installed successfully.

<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/mainGUI.png"  height="450" /></div>



## 2.Data inputting

The data would inputted in MRIES are three parts:
- Electrical stimulation-response data : Used to detect the CCEP response, should be *.edf* format.
- Reconstruction result : computed by Freesurfer.
- Subject infomation file : eg. age, sex, electrode implantation...
- Intracranial electrodes location file : computed by FIELD(引用).

This is the inputting file structure of the example data. 
<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/datainputall.png"  height="200" /></div>
Below we will take the example data to introduce these four parts in detail.

#### 1) Electrical stimulation-response data
Electrical stimulation-response data should be store in *edf* folder as *.edf* format. 
<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/edf.png"  height="200" /></div>

The electrode and channel information matched to *.edf* file should typed in the *edf_match.txt* file. 
<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/edf_match.png"  height="200" /></div>

#### 2) Reconstruction result

#### 3) Subject infomation file

#### 4) Intracranial electrodes location file

