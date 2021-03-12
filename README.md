Automate 16s
-------------------------------------------------
This project's aim is focused on automating processes of 16s rRNA analysis and learn the pipeline development tool Nextflow. The ultimate aim is to make the process more accessible. 

## Prerequisities
* Linux based system
* nextflow

## Install

```shell
$ git clone https://github.com/lorentzben/automate_16s_nf.git
```
### Note:
On setup, you will need to download one of the classifiers I generated below, or genertate your own. 
[515-806 Classifier](https://outlookuga-my.sharepoint.com/:u:/g/personal/bjl34716_uga_edu/EX3h7KrIg_5HqkUGhEkPyDYBbmEhBqsQzlLtIUAFyQzXDQ?e=Y2Y56B)
[Whole 16s Classifier](https://outlookuga-my.sharepoint.com/:u:/g/personal/bjl34716_uga_edu/EcydxUc2syZBnZLLT_Z9ISQBow7NdSRxYHxazvL9iCwxOQ?e=Yad1bX)
[Qiime Guide to Train Custom Classifer](https://docs.qiime2.org/2021.2/tutorials/feature-classifier/)
I have a script in this repo that will also generate the files above, however it takes about 12 hours to run. 

