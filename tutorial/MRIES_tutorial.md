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
Electrical stimulation-response data(*.edf* format) and *edf_match.txt* file should be stored in *edf* folder. 
<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/edf.png"  height="200" /></div>

The electrode and channel information matched to *.edf* file should typed in the *edf_match.txt* file. Every line corresponds to a *.edf* file. The fist row is the filename of *.edf* file. The second row is the stimulated electrode number of this file. The third row is the channel number corresponding to the electrode. 

For example, the "stimulation.edf" file contains response by channel 1-3 of electrode 7 and channel 2-6 of electrode 9. This line should be " stimulation [7,9] [1:3,2:6]". The picture below is the *edf_match.txt* file of example data. 
<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/edf_match.png"  height="100" /></div>

#### 2) Reconstruction result

The *label*，*mri*，*touch* and *surf* folders are reconstruction files by Freesurfer, and will be used in different visualization methods. Only copying these folders to the subject folder is OK.

#### 3) Subject infomation file
The *subjectinfo.txt* file should contain patient and electrode information. It should be the specfic format as the example data below.
<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/subjectinfo.png"  height="150" /></div>

- name: name of subject
- age : age of subject
- sex : sex of subject, female or male
- hemisphere : electrode implanted hemisphere, rh(right), lh(left) or bh(both).
- numelec: number of electrode.
- chan_array: channel number of every electrode.
- chan_name: label of every electrode.
#### 4) Intracranial electrodes location file

The *brain3D* folder is generated through the FIELD(引用). The file structure is shown below.


<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/brain3D.png"  height="200" /></div>

User can also use other electorde reconstruction toolboxs to obtain the location of electrode contacts and generate the *autocoordinate.mat* file without the FIELD by matlab. The structure of *autocoordinate.mat* is shown as below.


<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/autocoor.png"  height="200" /></div>

The *autocoordinate.mat* contain a matrix called *savecoor*. Every line of *savecoor* corresponds a contact. The first row is the number of contact. The secend row contain the electrode and contact information. eg. '107' means the 7th contact of first electrode, '1013' means the 13rd contact of 10th electrode. The 3rd-5th row is the *tkr* coordination of contact.



## 3.Processing pipline
After organizing data input, user can open the 'MRIES' main GUI to run the processing pipline. The processing pipline can be run by GUI button. Two mainly steps in this process,(1) subject and paramater setting and (2)  data processing.

<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/mainGUI.png"  height="450" /></div>

#### (1).Subject and paramater setting

Press the ‘Subject Information’ button, and subject's path eslection window will be shown. 

<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/subpath.png"  height="300" /></div>

Choose the subject folder ,then the ‘Electrode & Stimulation’ window will be shown. 

<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/elec_sti.png"  height="350" /></div>

The electrode label and channel number will be loaded. The stimulation setting should typed in this window, eg, number of stimulation pluse, bad contact, sampling rate and stimulation interval.

Press the 'Parameter Settings' button, and the parameter setting windown will be opened.User can set the parameters in this window. The specific meaning of the parameters can be seen in our future paper.

<div align=center><img src="https://github.com/SunKaijia0065/MRIES/blob/main/tutorial/image/elec_sti.png"  height="350" /></div>


