"""
author: zha
modified by xd
"""

import matplotlib.pyplot as plt 
import matplotlib.image as mpimg 
import numpy as np
from PIL import Image
import os
import csv

# path after manually labeling
reference_fig = r' figure roughlly matched previously '
# path before labeling where the ephys and probe figure pieced together
raw_fig = r' raw figure '

print('Click for user1\'s matching!')
probe0_after=Image.open(reference_fig)
plt.imshow(probe0_after)

user1_click_for_match=plt.ginput(30)  # order: ephys, probe; uppermostboundary, lowermostboundary, regionboundaries
print('user1\'s_click_for_match is %s' % user1_click_for_match)

length_user1_click_for_match=len(user1_click_for_match)

user1_y_ephys=[]
user1_y_histo_probe=[]

for j in range(length_user1_click_for_match):
    if j % 2==0:  # even for ephys 
        user1_y_ephys.append(user1_click_for_match[j][1])           
    else:   # odd for histology probe
        user1_y_histo_probe.append(user1_click_for_match[j][1])

print('user1_y_ephys is %s' % user1_y_ephys)
print('user1_y_histo_probe is %s' % user1_y_histo_probe)

#column_names_for_match=['user1_y_ephys','user1_y_histo_probe','user2_y_ephys','user2_y_histo_probe']
column_names_user1=['information','user1_y_ephys','user1_y_histo_probe']

row_content_user1=[os.path.splitext(raw_fig)[0],user1_y_ephys,user1_y_histo_probe]
print(row_content_user1)

(mousenumber,date,probenumber,point,AP,ML) = raw_fig.split('-')

savefilename_for_match= mousenumber + '-' +probenumber +'.csv'

savefilename_user1=os.path.splitext(raw_fig)[0] +'-user1.csv'

with open(savefilename_user1,'w',newline='')as csv_file:    # write csv,【user1's matching】
    filename_csv=csv.writer(csv_file)
    filename_csv.writerow(column_names_user1)
        
length_matchpoint_user1=len(row_content_user1[1])  # how many coordinates(clicks)
for k_user1 in range(length_matchpoint_user1):  # reorganize table
    row_content_for_match_list_user1=[row_content_user1[0],row_content_user1[1][k_user1],row_content_user1[2][k_user1]]
    with open(savefilename_user1,'a',newline='')as csv_file:     #write csv, mode 'a';【user1's matching】
        filename_csv=csv.writer(csv_file)
        filename_csv.writerow(row_content_for_match_list_user1)
print('user1 finished')

##################################################  user2's matching  ##################################################################

print('Click for user2\'s matching!')
probe0_after=Image.open(reference_fig)
plt.imshow(probe0_after)

user2_click_for_match=plt.ginput(30)   # order: ephys, probe; uppermostboundary, lowermostboundary, regionboundaries
print(user2_click_for_match)
length_user2_click_for_match=len(user2_click_for_match)

user2_y_ephys=[]
user2_y_histo_probe=[]

for j in range(length_user2_click_for_match):
    if j % 2==0:  # even for ephys
        user2_y_ephys.append(user2_click_for_match[j][1])
    else:   # odd for histology probe
        user2_y_histo_probe.append(user2_click_for_match[j][1])

print('user2_y_ephys is %s' % user2_y_ephys)
print('user2_y_histo_probe is %s' % user2_y_histo_probe)

column_names_user2=['user2_y_ephys','user2_y_histo_probe']
row_content_user2=[os.path.splitext(raw_fig)[0],user2_y_ephys,user2_y_histo_probe]
print(row_content_user2)
savefilename_user2=os.path.splitext(raw_fig)[0] +'-user2.csv'

with open(savefilename_user2,'w',newline='')as csv_file:    # write csv,【user2's matching】
    filename_csv=csv.writer(csv_file)
    filename_csv.writerow(column_names_user2)

length_matchpoint_user2=len(row_content_user2[1])  # how many coordinates(clicks)
for k_user2 in range(0,length_matchpoint_user2):  # reorganize table
    row_content_for_match_list_user2=[row_content_user2[0],row_content_user2[1][k_user2],row_content_user2[2][k_user2]]
    with open(savefilename_user2,'a',newline='')as csv_file:     #write csv, mode 'a';【user2's matching】
        filename_csv=csv.writer(csv_file)
        filename_csv.writerow(row_content_for_match_list_user2)
print('user2 finished')

####################################  calculate depth  ###############################################
if user1_y_ephys !=[]:
    user1_ephys_surf=user1_y_ephys[0]  
    user1_ephys_tip=user1_y_ephys[1]   
    user1_distance_ephys_3840=user1_ephys_tip-user1_ephys_surf # the real y-distance in the figure relative to 3840

    user1_y_ephys_3840=[]
    for each in user1_y_ephys:
        user1_y_ephys_3840.append((user1_ephys_tip-each)/user1_distance_ephys_3840*3840)
    print('user1_y_ephys_3840 is %s' % user1_y_ephys_3840)

if user2_y_ephys !=[]:
    user2_ephys_surf=user2_y_ephys[0]  
    user2_ephys_tip=user2_y_ephys[1]   
    user2_distance_ephys_3840=user2_ephys_tip-user1_ephys_surf # the real y-distance in the figure relative to 3840

    user2_y_ephys_3840=[]
    for each in user2_y_ephys:
        user2_y_ephys_3840.append((user2_ephys_tip-each)/user2_distance_ephys_3840*3840)
    print('user2_y_ephys_3840 is %s' % user2_y_ephys_3840)

##################################  combine results from users together  ###########################################

if user1_y_ephys ==[]:
    if user2_y_ephys==[]:
        print('user1 & user2 didn\'t do it !')  #  both users skipped
    else:  # user2 only
        column_names_for_match=['user2_y_ephys','user2_y_ephys_3840','user2_y_histo_probe']
        with open(savefilename_for_match,'w',newline='')as csv_file:    # table head 【user1 & user2's table】
                filename_csv=csv.writer(csv_file)
                filename_csv.writerow(column_names_for_match)
        for k_user2 in range(length_matchpoint_user2):
            row_content_for_match_list=[row_content_user2[1][k_user2],user2_y_ephys_3840[k_user2],row_content_user2[2][k_user2]] 
            with open(savefilename_for_match,'a',newline='')as csv_file:     # fill in the coordinates, mode 'a';【user1 & user2's table】
                filename_csv=csv.writer(csv_file)
                filename_csv.writerow(row_content_for_match_list)
        print('only user2 did !')

if user1_y_ephys != []:
    if user2_y_ephys==[]:  # user1 only
        column_names_for_match=['user1_y_ephys','user1_y_ephys_3840','user1_y_histo_probe']
        with open(savefilename_for_match,'w',newline='')as csv_file:    # table head 【user1 & user2's table】
                filename_csv=csv.writer(csv_file)
                filename_csv.writerow(column_names_for_match)
        for k_user1 in range(0,length_matchpoint_user1):
            row_content_for_match_list=[row_content_user1[1][k_user1],user1_y_ephys_3840[k_user1],row_content_user1[2][k_user1]] 
            with open(savefilename_for_match,'a',newline='')as csv_file:     # fill in the coordinates, mode 'a';【user1 & user2's table】
                filename_csv=csv.writer(csv_file)
                filename_csv.writerow(row_content_for_match_list)                
        print('only user1 did !')

    else: # both users did
        k_user1=length_matchpoint_user1
        k_user2=length_matchpoint_user2
        k_user_max = max(k_user1,k_user2)

        if k_user1<k_user2: # fill the empty coordinates with 0
            for k_addzero in range(k_user1,k_user2):
                user1_y_ephys.insert(k_addzero,0)
                user1_y_histo_probe.insert(k_addzero,0)
                user1_y_ephys_3840.insert(k_addzero,0)
        elif k_user1>k_user2: # fill the empty coordinates with 0
            for k_addzero in range(k_user2,k_user1):
                user2_y_ephys.insert(k_addzero,0)
                user2_y_histo_probe.insert(k_addzero,0)
                user2_y_ephys_3840.insert(k_addzero,0)
        else:
                column_names_for_match=['user1_y_ephys','user1_y_ephys_3840','user1_y_histo_probe','user2_y_ephys','user2_y_ephys_3840','user2_y_histo_probe']
                with open(savefilename_for_match,'w',newline='')as csv_file:    # table head 【user1 & user2's table】
                filename_csv=csv.writer(csv_file)
                filename_csv.writerow(column_names_for_match)
                for k in range(k_user_max):
                row_content_for_match_list=[row_content_user1[1][k],user1_y_ephys_3840[k],row_content_user1[2][k],row_content_user2[1][k],user2_y_ephys_3840[k],row_content_user2[2][k]] 
                with open(savefilename_for_match,'a',newline='')as csv_file:     # fill in the coordinates, mode 'a';【user1 & user2's table】
                        filename_csv=csv.writer(csv_file)
                        filename_csv.writerow(row_content_for_match_list)            
print( 'match finished')
            






