# Welcome to the MRIES tutorial
This tutorial is an introduction to [MRIES](https://github.com/SunKaijia0065/MRIES_tutorial) processing and visualization.



The MRIES Tutorial is split into five parts, [installing MRIES](#1)
, [data inputting](#2), 
[processing pipline](#3),
 [visualization](#4) 
and [result outputting](#5).
We will take the [example data](https://figshare.com/articles/dataset/MRIES_s_sample_data/14376110) as example to show the use of MRIES.



<div id="1"></span>

## 1. Installing MRIES: 
1. Start Matlab.

2. Add the *MRIES* folder to path of Matlab. 

<img  src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial/main/image/path.png"/>

3. Type "MRIES" at the Matlab command prompt ">>" and press enter. If MRIES main GUI will be shown, MRIES will be installed successfully.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/mainGUI.png"  height="450" /></div>






<div id="2"></span>

## 2.Data inputting

The data would inputted in MRIES are three parts:
- Electrical stimulation-response data : Used to detect the CCEP response, should be *.edf* format.
- Reconstruction result : computed by Freesurfer.
- Subject infomation file : eg. age, sex, electrode implantation...
- Intracranial electrodes location file : computed by FIELD[(Qin 2017)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5314105/pdf/fninf-11-00010.pdf).

This is the inputting file structure of the example data. 
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/datainputall.png"  height="200" /></div>
Below we will take the example data to introduce these four parts in detail.

#### 1) Electrical stimulation-response data
Electrical stimulation-response data(*.edf* format) and *edf_match.txt* file should be stored in *edf* folder. 
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/edf.png"  height="200" /></div>

The electrode and channel information matched to *.edf* file should typed in the *edf_match.txt* file. 
Every line corresponds to a *.edf* file. 
The fist row is the filename of *.edf* file. The second row is the stimulated electrode number of this file. The third row is the channel number corresponding to the electrode. 

For example, the "stimulation.edf" file contains response by channel 1-3 of electrode 7 and channel 2-6 of electrode 9. This line should be " stimulation [7,9] [1:3,2:6]". The picture below is the *edf_match.txt* file of example data.
 
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/edf_match.png"  height="100" /></div>

#### 2) Reconstruction result

The *label*，*mri*，*touch* and *surf* folders are reconstruction files by Freesurfer, and will be used in different visualization methods. Only copying these folders to the subject folder is OK.

#### 3) Subject infomation file
The *subjectinfo.txt* file should contain patient and electrode information. It should be the specfic format as the example data below.
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/subjectinfo.png"  height="150" /></div>

- name: name of subject
- age : age of subject
- sex : sex of subject, female or male
- hemisphere : electrode implanted hemisphere, rh(right), lh(left) or bh(both).
- numelec: number of electrode.
- chan_array: channel number of every electrode.
- chan_name: label of every electrode.
#### 4) Intracranial electrodes location file

The *brain3D* folder is generated through the FIELD(引用). The file structure is shown below.


<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/brain3D.png"  height="200" /></div>

User can also use other electorde reconstruction toolboxs to obtain the location of electrode contacts and generate the *autocoordinate.mat* file without the FIELD by matlab. The structure of *autocoordinate.mat* is shown as below.


<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/autocoor.png"  height="200" /></div>

The *autocoordinate.mat* contain a matrix called *savecoor*. Every line of *savecoor* corresponds a contact. The first row is the number of contact. The secend row contain the electrode and contact information. eg. '107' means the 7th contact of first electrode, '1013' means the 13rd contact of 10th electrode. The 3rd-5th row is the *tkr* coordination of contact.


<div id="3"></span>

## 3.Processing pipline
After organizing data input, user can open the 'MRIES' main GUI to run the processing pipline. The processing pipline can be run by GUI button. Two mainly steps in this process,(1) subject and paramater setting and (2)  data processing.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/mainGUI.png"  height="450" /></div>

#### (1).Subject and paramater setting

Press the ‘Subject Information’ button, and subject's path eslection window will be shown. 

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/subpath.png"  height="300" /></div>

Choose the subject folder ,then the ‘Electrode & Stimulation’ window will be shown. 

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/elec_sti.png"  height="350" /></div>

The electrode label and channel number will be loaded. The stimulation setting should typed in this window, eg, number of stimulation pluse, bad contact, sampling rate and stimulation interval.

Press the 'Parameter Settings' button, and the parameter setting windown will be opened as below.User can set the parameters in this window. The specific meaning of the parameters can be seen in our future paper.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/setting.png"  height="450" /></div>

Subject information and parameter setting can be shown in command line.
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/setting_output.png"  height="300" /></div>

#### (2).Data processing

Data processing can be run as complete or separate model.

Pressins the 'Data Processing pipeline' button can dirctly run whole pipeline of processing.
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/pipelinebutton.png"  width="300" /></div>

User can also press the buttons in ‘Separate Propressing Function’ panel to run the processing separatly. 
The detailed explanation of processing steps can be seen in our future paper.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/sparatebutton.png"  width="300" /></div>

Processing detail can be shown in command line.
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/lineall.png"  width ="400" /></div>

- Tips: There are two case in epoching step, With and without marker channel. When the data contain the marker channel, user only need type label of the marker channel in parameter setting window. When the data don't contain the marker channel, user neeed empty the marker channel and the command will remind the user to type the channel number used to detect in the processing.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/lineepoch.png"  width ="400" /></div>



<div id="4"></span>

## 4. Visualization
Press the 'Visulization' button can open the MRIES's visulization GUI, MRIESviewer as below.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/MRIESviewer.png"  height="400" /></div>

The subject and electrode information can be shown in the .
The connectivity method and stimulation contact to show can be chosen in the .
The *autocoordinate.mat* file should be input by the 'coordinates' button.

There are four visualization methods in MRIES.
(1)Matrix, (2)Circle map, (3)Volume and (4)Surface. 
After choosing the method and stimilation to display, user can have the four kinds of visualization as below. 
We use the example data to display these methods. 
The detailed explanation of visualization can be seen in our future paper.


#### (1) Matrix
The controller ‘Conn Mat’ in visualization GUI can open the matrix GUI. The connectivity can be shown as hot-map. The threshold can change by the Slider.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/Matrix.png"  height="400" /></div>

#### (2) Circle map

The controller ‘Circle Map’ in visualization GUI can open the circle map GUI.
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/circle.png"  height="400" /></div>

#### (3) Volume

The volume visualization can be shown in MRIESviewer GUI. This method merge the MR image, CT image and response result.
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/volume.png"  height="400" /></div>

#### (4) Surface

The controller ‘Surface’ in visualization GUI can open the Surface GUI. User can rotate the view in this GUI. 

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/surface.png"  height="400" /></div>


<div id="5"></span>

## 5. Result outputting


After all processsing, the *data*, *stimulation* and *Matrix* folders will be generated  in the subject folder by MRIES.
<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/folderresult.png"  height="200" /></div>

The *data* folder contains raw signal separated into the different stimulated contact-pairs. For example, '*ccep_elec_14_15.mat*' means the data by stimulating 14th and 15th channel.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/folderdata.png"  height="200" /></div>

Every file in *stimulationdata* folder includes all response indicators from one stimulation contact-pairs. 
For example, '*ccep_elec_B2_B3_All.mat*' means the response detection result by stimulating 2nd and 3rd contact of 'B' electrode.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/folderstimu.png"  height="200" /></div>

The *Matrix* folder includes all the connectivity matrix of all kinds of indicators, which can be used in the visualization GUIs. 
For example, '*conn_matrix_LF_RMS.mat*' means the low-frecquency RMS matrix. These matrix can be used to further analysis and calculation.

<div align=center><img src="https://raw.githubusercontent.com/SunKaijia0065/MRIES_tutorial//main/image/folderMatrix.png"  height="200" /></div>



























