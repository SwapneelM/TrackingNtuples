#!/usr/bin/env python
# coding: utf-8

# # Table of Contents
# 
# 
# ### Part I - Preprocessing
# ---------------------------
# 
# 1. [View Keys in Root Data](#View-the-Keys-in-the-Imported-Data)
# 
# 2. [Optimisation Tests I: Python](#Optimisation-Tests)
# 
# 3. [Data Conversion using Uproot](#Load-Data-into-Arrays)
# 
# 4. Preprocessing
# 	- [Optimisation](#Preprocessing-1:-Reformat-List-of-Indices-to-Sets-of-Indices-for-each-Rechit)
# 	- [Convert to Dataframes](#Preprocessing-2:-Add-all-data-into-dataframes)
# 	- [Create Global Rechit Dataframe](#Create-a-Global-Dataframe-of-Rechits)
# 
# 5. [Rechit to Track Matching](#Match-the-Rechits-to-Tracks-and-Create-a-Global-Array-of-Tracks)
# 
# 6. [Optimisation Tests II: Dataframes](#Test-Performance-of-df.loc-versus-multi-index-retrieval-[-][-])
# 
# 
# ### Part II: Raw Data Analysis and Plots
# -----------------------------------------
# 
# 1. [Count Data in Track to Rechit Map](#Analyse-the-data-stored-in-the-track_to_rechit_map_)
# 
# 2. [Generate Match Count Plots](#Generate-Plots)
# 	- [Counting Matched vs. Unmatched Rechits](#Analyse-Matched/Unmatched-Rechits)
# 	- [Plot Count of Track and TP Matched Rechits](#Plot-the-Rechits-Matched-to-TP,-Track,-or-Unmatched)
# 
# 3. Track Analysis
# 	- [Compare Eta between Tracks and their Rechits](#Compare-Track-and-Matched-Rechit-Eta)
# 	- [Plot Track Parameter Distribution Histograms](#Plot-Track-Parameters)
# 	- [Analyse/Filter High-Pt Events](#Filter-High-Pt-Events)
# 
# 4. Rechit/Simhit Analysis
# 	- [Plot Simhit Distribution 2D](#Plot-SimHit-Distribution-in-X-and-Y-Axes)
# 	- [Simhit Match Count](#Count-Simhits-Matched-to-Tracks)
# 	- [Realistic Geometry Simulation: Hole](#Plot-the-Hole-in-the-Data-(2D-Plot;-3D-Axes))
# 	- [Plot MonoRechit Distribution 2D](#Visualize-the-Mono-Rechits)
# 	- [Plot Rechit Parameter Distribution Histograms](#Plot-Rechit-Parameters)
# 	- [Plot StereoRechit Distribution 2D](#Visualize-the-Stereo-Rechits)
# 
# 5. Data Storage
# 	- [Verify Data is not Corrupted](#Testing-Integrity-of-internal-data-storage)
# 	- [Write Data to Serialized Output Format](#Data-Storage-for-TF/PyTorch/Graph-Library)
# 
# 
# ### Part III: Filtered (Cut) Data Analysis and Plots
# --------------------------------------------------
# 
# 1. [Place Cuts on Rechits](#Place-the-Cuts-on-Rechits-=>-eta-(-0.9,-0.9))
#     - [Scatter Plot of Filtered Rechits](#Plot-the-Hits-without-Connections)
# 
# 2. [Plot 2D Rechit Parameters](#Plot-the-2D-Rechit-Parameters-for-Filtered-Hits)
# 
# 3. [Place Cuts on Tracks](#Place-the-cuts-on-tracks-by-Eta-and-Pt)
# 
# 4. [Plot only Filtered Rechits](#Plot-only-the-filtered-rechits)
# 
# 5. [Plot Matched/Unmatched Track Distribution](#Plot-Track-Distribution)

# In[4]:


get_ipython().magic(u'matplotlib notebook')
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import uproot
import pandas as pd
from collections import OrderedDict


# ## View the Keys in the Imported Data

# In[5]:


gen_event_ = "ttbar-100"
number_of_events_ = 100
outfile_ = "outfile-" + gen_event_ + ".root"
data_ = uproot.open(outfile_)["ntuples"]["tree"]
# data_.keys()


# ## Check the Integrity of the Imported Data 

# In[6]:


stereo_tp_idx_ = data_.array('stereoTPIndex')
mono_tp_idx_ = data_.array('monoTPIndex')
track_tp_idx_ = data_.array('trackTPIndex')

# Check that both have been generated for the same number of events
# Just for clarity
assert len(track_tp_idx_) == len(stereo_tp_idx_), "Track and Stereo Number of Events do not match"
assert len(track_tp_idx_) == len(mono_tp_idx_), "Track and Mono Number of Events do not match"
print ("\nTotal", len(track_tp_idx_), "events")


# # Optimisation Tests

# In[7]:


def list_to_set(input_array_):
    '''
    Format: 3-level nested lists - [[[...] ...] ...]
    '''
    output_array_ = []
    for index_ in range(len(input_array_)):
        output_array_.append([])
        for second_list_ in input_array_[index_]:
            output_array_[index_].append(set(second_list_))
    return output_array_


# In[8]:


mono_tp_idx_set_ = list_to_set(mono_tp_idx_)


# ## Load Data into Arrays

# In[9]:


'''
Load the track parameters into the respective arrays to be added into the rechit_param_global dataframe
'''

rechit_cartesian_ = OrderedDict({})
for key in ['stereoHitX', 'stereoHitY', 'stereoHitZ', 'monoHitX', 'monoHitY', 'monoHitZ']:
    rechit_cartesian_[key] = data_.array(key)

rechit_polar_ = OrderedDict({})
for key in ['stereoHitR', 'stereoHitEta', 'stereoHitPhi', 'monoHitR', 'monoHitEta', 'monoHitPhi']:
    rechit_polar_[key] = data_.array(key)


# ## Preprocessing 1: Reformat List of Indices to Sets of Indices for each Rechit

# In[10]:


# Convert all tracking particle index lists to sets for faster search

mono_tp_idx_set_ = list_to_set(mono_tp_idx_)
stereo_tp_idx_set_ = list_to_set(stereo_tp_idx_)
track_tp_idx_set_ = list_to_set(track_tp_idx_)


# ## Preprocessing 2: Add all data into dataframes

# In[11]:


import pandas as pd
from pandas import DataFrame as df


# ### Create a Global Dataframe of Rechits

# In[12]:


'''
Adding stereo and mono rechit data into a global dataframe

:event_id: int
:rechit_id: int
:track_id: int
:rechit_ids: list(int)
:track_ids: list(int)
:track_tp_index: set(int)  # iterating over sets has lower complexity
:rechit_tp_index: set(int)  # iterating over sets has lower complexity
:match_count: int  # count the number of rechits/tracks matched to the given track/rechit
:rechit_tp_index_: event-based list of rechit-based list of sets of int (tp_index)

'''
def create_global_rechit_df(stereo_tp_idx_, mono_tp_idx_, rechit_cartesian_dict_, rechit_polar_dict_):
    rechit_global_map_ = OrderedDict({'event_id': [], 'rechit_id': [], 'rechit_tp_index': [],
                                      'track_ids': [], 'match_count': [], 'rechit_local_id': []})
    rechit_param_global_map_ = OrderedDict({'event_id': [], 'rechit_id':[], 'rechit_x': [], 'rechit_y': [], 'rechit_z': [], 
                                            'rechit_r': [], 'rechit_phi': [], 'rechit_eta': [], 'rechit_local_id': []})
    global_counter_ = 0
    
    if len(stereo_tp_idx_) != len(stereo_tp_idx_):
        raise ValueError('Rechit arrays represent differing event lengths [stereo, mono]:', len(stereo_tp_idx_), len(mono_tp_idx_))
    
    for event_id_ in range(len(stereo_tp_idx_)):
        # Count the number of rechits in that event
        event_rechit_count_ = len(stereo_tp_idx_[event_id_]) + len(mono_tp_idx_[event_id_])

        rechit_global_map_['event_id'].extend([event_id_] * event_rechit_count_)  
        # appends SAME instance of [event_id] event_rechit_count_ times
        
        rechit_global_map_['rechit_id'].extend(
            range(global_counter_, global_counter_ + event_rechit_count_))     
        rechit_global_map_['rechit_tp_index'].extend(stereo_tp_idx_[event_id_])
        rechit_global_map_['rechit_tp_index'].extend(mono_tp_idx_[event_id_])
        rechit_global_map_['track_ids'].extend([[] for _ in range(event_rechit_count_)])
        rechit_global_map_['match_count'].extend([0 for _ in range(event_rechit_count_)])
        rechit_global_map_['rechit_local_id'].extend(range(event_rechit_count_))
        
        # Extend the hit_param_global_map_ with rechit parameters
        rechit_param_global_map_['rechit_id'].extend(
            range(global_counter_, global_counter_ + event_rechit_count_))
        rechit_param_global_map_['event_id'].extend([event_id_] * event_rechit_count_)  
        rechit_param_global_map_['rechit_x'].extend(rechit_cartesian_dict_['stereoHitX'][event_id_])
        rechit_param_global_map_['rechit_x'].extend(rechit_cartesian_dict_['monoHitX'][event_id_])
        rechit_param_global_map_['rechit_y'].extend(rechit_cartesian_dict_['stereoHitY'][event_id_])
        rechit_param_global_map_['rechit_y'].extend(rechit_cartesian_dict_['monoHitY'][event_id_])
        rechit_param_global_map_['rechit_z'].extend(rechit_cartesian_dict_['stereoHitZ'][event_id_])
        rechit_param_global_map_['rechit_z'].extend(rechit_cartesian_dict_['monoHitZ'][event_id_])
        
        rechit_param_global_map_['rechit_r'].extend(rechit_polar_dict_['stereoHitR'][event_id_])
        rechit_param_global_map_['rechit_r'].extend(rechit_polar_dict_['monoHitR'][event_id_])
        rechit_param_global_map_['rechit_phi'].extend(rechit_polar_dict_['stereoHitPhi'][event_id_])
        rechit_param_global_map_['rechit_phi'].extend(rechit_polar_dict_['monoHitPhi'][event_id_])
        rechit_param_global_map_['rechit_eta'].extend(rechit_polar_dict_['stereoHitEta'][event_id_])
        rechit_param_global_map_['rechit_eta'].extend(rechit_polar_dict_['monoHitEta'][event_id_])
        rechit_param_global_map_['rechit_local_id'].extend(range(event_rechit_count_))
        global_counter_ += event_rechit_count_
    # Convert dict to dataframe
    rechit_global_df_ = df.from_dict(rechit_global_map_)
    rechit_param_global_df_ = df.from_dict(rechit_param_global_map_)
    return rechit_global_df_, rechit_param_global_df_
    
# Check Memory Usage of DataFrame
# print rechit_global_df_.memory_usage(deep=True)
# print rechit_param_global_df_.memory_usage(deep=True)


# In[14]:


'''
Create the Global Rechit Array and Global Rechit Parameters Array'''
rechit_global_df_uncut_, rechit_param_global_df_uncut_ = create_global_rechit_df(
    stereo_tp_idx_, mono_tp_idx_, rechit_cartesian_, rechit_polar_)
#print rechit_global_df_.head(10)


# ## Place the Cuts (create DF for Graph Networks)

# In[19]:


'''Check the maximum number of hits in an event'''
max_len_ = 0 
for i in range(number_of_events_):
    len_idx_ = len(stereo_tp_idx_[i]) + len(mono_tp_idx_[i])
    if len_idx_ > max_len_:
        max_len_ = len_idx_
        #print max_len_
        
print ("Maximum hits in an event are: ", max_len_)
     


# In[20]:


'''This is done here to generate a reduced number of local indices for the tracks to match to rechits.
We will replace the rechit_global_df_ generated above and used below for matches with this new dataframe.'''

intermediate_df_ = rechit_param_global_df_uncut_[rechit_param_global_df_uncut_['rechit_eta'] <= 0.9]
rechit_param_global_df_ = intermediate_df_[intermediate_df_['rechit_eta'] >= -0.9].copy()
rechit_global_df_ = rechit_global_df_uncut_.iloc[rechit_param_global_df_['rechit_id']].copy()

# Reset the Index of the Cut Dataframe that will become the new Global DataFrame
# This will lose the former global rechit index - can this affect the analysis in the future?
rechit_global_df_.index = pd.RangeIndex(len(rechit_global_df_.index))  
rechit_param_global_df_.index = pd.RangeIndex(len(rechit_global_df_.index))  

# Reset the local_rechit_ids for graph networks to have sequential nodes
# And so that the node feature vector can be simpler to create sequentially
rechit_local_id_dict_ = {'rechit_local_id' : []}
# Find the minimum number of rechits in the final list of events
min_num_of_rechits_ = 9999
for event_id_ in range(number_of_events_):
    # Retrieve the subset of the global rechit dataframe for this event_id
    rechit_local_range_ = range(len(rechit_global_df_[rechit_global_df_['event_id']==event_id_]))
    rechit_local_id_dict_['rechit_local_id'].extend(rechit_local_range_)
    if rechit_local_range_[-1] < min_num_of_rechits_:
        min_num_of_rechits_ = rechit_local_range_[-1]

# Update the Global Rechit IDs
rechit_global_id_dict_ = {}
rechit_global_id_dict_['rechit_id'] = range(len(rechit_global_df_))
rechit_global_df_.update(pd.DataFrame.from_dict(rechit_local_id_dict_))    
rechit_param_global_df_.update(pd.DataFrame.from_dict(rechit_local_id_dict_))    

# Update the Local Rechit IDs
rechit_global_df_.update(pd.DataFrame.from_dict(rechit_global_id_dict_))    
rechit_param_global_df_.update(pd.DataFrame.from_dict(rechit_global_id_dict_))

print (len(rechit_param_global_df_), "of", len(rechit_param_global_df_uncut_), float(len(rechit_param_global_df_))/float(len(rechit_global_df_uncut_)), "hits remain")


# In[ ]:





# ## Match the Rechits to Tracks and Create a Global Array of Tracks

# In[23]:


'''
Match Rechits to Tracks.
Create the Global Track Array and Global Track Parameter Array.
'''
# TODO: Refactor this to enable placing track cuts before forming dataframe and reduce processing by 75%
# The 75% metric follows from: For 100 events track cuts reduce tracks by 75%
# Initialize the Global Track Parameter Map
track_param_global_map_ = OrderedDict({})
for key in ['track_id', 'track_eta', 'track_phi', 'track_qoverp', 'track_dxy', 'track_dsz', 'track_pt']:
    track_param_global_map_[key] = []
    
# Define the dictionaries to be cast into dataframes
track_to_rechit_map_ = OrderedDict({'event_id': [], 'track_id': [], 'track_tp_index': [], 
                                    'rechit_ids': [], 'match_count': [], 'rechit_local_ids': []})

# Future Requirement?
rechit_to_track_map_ = OrderedDict({'event_id': [], 'rechit_id': [], 'rechit_tp_index': [],
                                    'track_ids': [], 'match_count': []})

# Initialize the Global Track ID
global_track_id_ = 0

for event_id_ in range(len(track_tp_idx_)):
    
    num_tracks_in_event_ = len(track_tp_idx_[event_id_])

    # Add track data to the dict in an efficient manner
    track_to_rechit_map_['event_id'].extend([event_id_] * num_tracks_in_event_)
    
    global_track_id_range_ = range(global_track_id_, global_track_id_ + num_tracks_in_event_)
    
    track_to_rechit_map_['track_id'].extend(global_track_id_range_)
    track_to_rechit_map_['track_tp_index'].extend(track_tp_idx_[event_id_])
    
    # Append multiple empty lists in place of the values not filled yet
    track_to_rechit_map_['match_count'].extend([] for _ in range(num_tracks_in_event_))
    track_to_rechit_map_['rechit_ids'].extend([] for _ in range(num_tracks_in_event_))
    track_to_rechit_map_['rechit_local_ids'].extend([] for _ in range(num_tracks_in_event_))
    
    # Fill in the Global Track Parameters
    track_param_global_map_['track_id'].extend(global_track_id_range_)
    track_param_global_map_['track_eta'].extend(data_.array('trackEta')[event_id_])
    track_param_global_map_['track_phi'].extend(data_.array('trackPhi')[event_id_])
    track_param_global_map_['track_pt'].extend(data_.array('trackPt')[event_id_])
    track_param_global_map_['track_qoverp'].extend(data_.array('qoverp')[event_id_])
    track_param_global_map_['track_dxy'].extend(data_.array('dxy')[event_id_])
    track_param_global_map_['track_dsz'].extend(data_.array('dsz')[event_id_])
    
    # Retrieve the subset of the global rechit dataframe for this event_id
    event_df_ = rechit_global_df_[rechit_global_df_['event_id']==event_id_]
    
    # Check the TPs matched to tracks and find rechits for each TP (Stereo and Mono)
    for track_tp_list_ in track_tp_idx_[event_id_]:
        rechit_matches_ = []
        rechit_local_matches_ = []
        if len(track_tp_list_) == 0:
            continue
            
        if len(track_tp_list_) == 1:

            # Iterate over the index and values of each rechit tp index list
            for (idx_, tp_idx_list_) in event_df_['rechit_tp_index'].iteritems():
                # Find the match for the first tp index in the track tp list
                if track_tp_list_[0] in tp_idx_list_:
                    rechit_matches_.append(event_df_.loc[idx_, 'rechit_id'])
                    rechit_local_matches_.append(event_df_.loc[idx_, 'rechit_local_id'])
                    # Append the global track id to the rechit
                    event_df_.loc[idx_, 'track_ids'].append(global_track_id_)
            track_to_rechit_map_['match_count'][global_track_id_] = len(rechit_matches_)
            track_to_rechit_map_['rechit_ids'][global_track_id_] = set(rechit_matches_)
            track_to_rechit_map_['rechit_local_ids'][global_track_id_] = set(rechit_local_matches_)
            
        # If track has multiple tp indices, pick the one with the most hits

        # Note: This approach *possibly* creates match issues if the tp index with more rechit matches
        # has more 'common' hits with other tracks and is later discarded due to the common hits 
        # belonging to other tracks
        if len(track_tp_list_) > 1:
            rechit_matches_array_ = []
            rechit_local_matches_array_ = []
            match_count_array_ = []
            
            print ("Found multiple TP indices in event", event_id_, "for global track") 
            print (global_track_id_, track_tp_list_)
            
            for track_idx_ in track_tp_list_:
                rechit_matches_ = []
                rechit_local_matches_ = []
                
                # Iterate over the index and values of each rechit tp index list
                for (idx_, tp_idx_list_) in event_df_['rechit_tp_index'].iteritems():
                    if track_idx_ in tp_idx_list_:
                        rechit_matches_.append(event_df_.loc[idx_,'rechit_id'])
                        rechit_local_matches_.append(event_df_.loc[idx_,'rechit_local_id'])
                        # Append the global track id to the rechit
                        event_df_.loc[idx_, 'track_ids'].append(global_track_id_)
                rechit_matches_array_.append(rechit_matches_)
                rechit_local_matches_array_.append(rechit_local_matches_)
                match_count_array_.append(len(rechit_matches_))
            
            # Store the global rechit ids and count of matches in a temporary list
            for key, value in zip(match_count_array_, rechit_matches_array_):
                tmp_dict_.append((key, value))
            
            # Pick the largest number of matches and corresponding global rechit ids
            tmp_dict_ = sorted(tmp_dict_, reverse=True)
            track_to_rechit_map_['match_count'][global_track_id_] = tmp_dict_[0][0]
            track_to_rechit_map_['rechit_ids'][global_track_id_] = tmp_dict_[0][1]
            track_to_rechit_map_['rechit_local_ids'][global_track_id_] = set(rechit_local_matches_array_)

        # Check duplicates
        if len(set(rechit_matches_)) < len(rechit_matches_):
            raise ValueError('rechit_matches_ has duplicate values: Some Rechits are being matched twice!')
        
        # Increment the Global Track ID
        global_track_id_ += 1
    rechit_global_df_.update(event_df_, join='left')
    track_param_global_df_ = df.from_dict(track_param_global_map_)
track_global_df_ = df.from_dict(track_to_rechit_map_)

#Update the match_count for rechits based on the number of total matched tracks
match_count_tmp_dict_ = OrderedDict({'match_count': [len(track_id_list_) for track_id_list_ in rechit_global_df_['track_ids']]})
print ("Maximum tracks matched for one particle:", max(match_count_tmp_dict_['match_count']))

rechit_global_df_.update(df.from_dict(match_count_tmp_dict_))


# ## Analyse the data stored in the track_to_rechit_map_

# In[25]:


track_to_rechit_df_ = df.from_dict(track_to_rechit_map_)
#print track_to_rechit_df_[track_to_rechit_df_['event_id']==11].head(10)

# Calculate the average number of hits per track
average_rechits_per_track_ = 0
len_array_ = []
for rechit_list_ in track_to_rechit_df_['rechit_ids']:
    average_rechits_per_track_ += len(rechit_list_)
    len_array_.append(len(rechit_list_))

print ("Average Rechits per track:", average_rechits_per_track_/len(track_to_rechit_df_['rechit_ids']))
print ("Max. matched hits to track:", max(len_array_), "; Global track id:", len_array_.index(max(len_array_)))


# Test to check if the correct tp index has been matched
# Change the value of 'trk_id_' to any track that you know has some hits
trk_id_ = len_array_.index(max(len_array_))
# print track_to_rechit_df_.loc[trk_id_]
for rechit_id in track_to_rechit_df_.loc[trk_id_]['rechit_ids']:
    for track_idx_ in track_to_rechit_df_.loc[trk_id_]['track_tp_index']:
        if track_idx_ in rechit_global_df_.loc[rechit_id]['rechit_tp_index']:
            continue
        else:
            print ("Error: Track and rechit TP index does not match!")
            break


# # Generate Plots

# In[26]:


from mpl_toolkits.mplot3d import Axes3D
from cycler import cycler
from matplotlib.colors import Colormap

#fig_ = plt.figure()
#ax_ = Axes3D(fig_)


# ### Analyse Matched/Unmatched Rechits

# In[27]:


'''
Count the number of matched, unmatched, and total rechits/tracks in the dataframe (PER EVENT)

Store the count of unmatched, tp_matched, track/rechit_matched, and total rechits/tracks PER EVENT in an array of length number_of_events_
Store all four above arrays (unmatched, tp_matched, track/rechit_matched, total) in a dictionary
'''

def count_matched_items(item_type_):
    other_item_ids_ = 'track_ids' if (item_type_=='rechit') else 'rechit_ids'
    other_item_matched_ = 'track_matched' if (item_type_=='rechit') else 'rechit_matched'
    item_id_ = item_type_ + '_id'
    item_tp_index_ = item_type_ + '_tp_index'

    # Initialize one array for counts and one for ids of matched/unmatched rechits
    item_count_dict_ = OrderedDict({other_item_matched_:[], 'unmatched':[], 'tp_matched':[], 'total':[]})
    item_id_dict_ = OrderedDict({'tp_matched':[], other_item_matched_:[], 'unmatched':[]})

    for event_id_ in range(number_of_events_):
        
        # Create a slice of the dataframe with the data for that event
        event_df_ = (rechit_global_df_[rechit_global_df_['event_id']==event_id_]) if (item_type_=='rechit') else (track_global_df_[track_global_df_['event_id']==event_id_])

        # Count the number of matched, unmatched, and total rechits 
        num_matched_ = sum(event_df_['match_count'] > 0)
        num_unmatched_ = sum(event_df_['match_count'] == 0)
        num_total_ = event_df_.shape[0]  # number of rows/rechits in the event
        
        # Find and store the indices of matched and unmatched rechits
        
        item_id_dict_[other_item_matched_].append(set(event_df_.loc[event_df_['match_count'] > 0, (item_id_)].tolist()))
        item_id_dict_['unmatched'].append(set(event_df_.loc[event_df_['match_count'] == 0, (item_id_)].tolist()))
        
        # Sanity checks to ensure data has been added into the dataframe corrrectly
        assert num_total_ == (num_matched_ + num_unmatched_),         "Rechit counts (unmatched, matched, total) do not add up"
            
        if item_type_ == 'rechit':
            # Check the number of total rechits for the event is the same as in raw data
            assert (len(rechit_global_df_[rechit_global_df_['event_id']==event_id_])) == num_total_,             "Rechits in dataframe %d and stereo_tp_idx_ %d do not match" % (num_total_, len(stereo_tp_idx_[event_id_]))
        
        elif item_type_ == 'track':
            # Check the number of total tracks for the event is the same as in raw data
            assert len(track_tp_idx_[event_id_]) == num_total_,             "Tracks in dataframe %d and track_tp_idx_ %d do not match" % (num_total_, len(track_tp_idx_[event_id_]))
    
        # Append the hit counts into the dataframe
        item_count_dict_['unmatched'].append(num_unmatched_)
        item_count_dict_['total'].append(num_total_)
        
        # TODO: Why is default value for tracks -2 and rechits None?
        # Criteria for tracks is to check if -2 is in the track_tp_index
        # Because default match to tp index value is -2
        if item_type_ == 'track':
            tp_criteria_ = [(-2 not in list_) for list_ in event_df_[item_tp_index_]]
        
        # Criteria for rechits is to check if length of rechit_tp_index is greater than 0
        # Because default match to tp index is none
        elif item_type_ == 'rechit':
            tp_criteria_ = [(len(list_) > 0) for list_ in event_df_[item_tp_index_]]
            #print len(event_df_[tp_criteria_])
        
        item_count_dict_['tp_matched'].append(len(event_df_[tp_criteria_]))
        item_id_dict_['tp_matched'].append(event_df_[tp_criteria_])
        
        # Criteria for filtering rechits matched to tracks based on 'track_ids' column
        other_item_criteria_ = [len(list_) > 0 for list_ in event_df_[other_item_ids_]]
        item_count_dict_[other_item_matched_].append(len(event_df_[other_item_criteria_]))
    
    return item_count_dict_, item_id_dict_

def plot_matched_vs_unmatched(item_count_, keys_, item_type_):
    ax_ = plt.subplot()
    alpha_ = 0.4
    for key in keys_:
        ax_.hist(item_count_[key], histtype='stepfilled', bins=number_of_events_, 
             orientation='vertical', alpha=alpha_, label=key)
        alpha_ += 0.2
    plt.grid(True)
    plt.legend()
    plt.ylabel('Frequency')
    plt.xlabel('Count of ' + item_type_)
    plt.title(item_type_ + ' Distribution')
    plt.savefig('plots/' + gen_event_ + '/' + item_type_ + '/matchdistribution')
    plt.show()
    return


# In[28]:


track_count_, track_ids_ = count_matched_items('track')
rechit_count_, rechit_ids_ = count_matched_items('rechit')


# ### Testing Integrity of internal data storage 

# In[30]:


# Correlate the data to confirm the dataframe has not been corrupted
hit_tp_count_ = {}

for (id_, tp_idx_list_) in rechit_global_df_["rechit_tp_index"].iteritems():
    tp_len_ = len(tp_idx_list_)
    if tp_len_ in hit_tp_count_:
        hit_tp_count_[tp_len_] += 1
    else:
        hit_tp_count_[tp_len_] = 1
print (hit_tp_count_)


# ## Data Storage for TF-DeepHGCal/PyTorch/Graph Library

# In[31]:


'''DataFrame Documentation for Pandas states that writing and reading from msgpack is an experimental feature.
It is to be released soon, but please use it with care to ensure data is not corrupted.

Note: When working with large datasets (>1000 events), you will not be able to save the data.
The filesize for rechit_global_df_ is 76 MB for 100 events thus 760 MB for 1000 events and so on.'''


# In[32]:


'''
Writing to serialized format fails in case of copied dataframes as the columns are sets

Solution: Iterate over all Global DataFrames, find the columns to replace, 
and replace with lists instead of sets so that they are serializable
'''
'''
for dataframe_ in [track_global_df_, track_param_global_df_, rechit_global_df_, rechit_param_global_df_]:
    dataframe_to_update_ = dataframe_.copy(deep=True)
    columns_to_replace_ = ['rechit_ids', 'rechit_local_ids', 'rechit_tp_index', 'track_tp_index', 'track_matches']
    for column_name_ in columns_to_replace_:    
        if column_name_ in dataframe_to_update_:
            list_arr_ = []
            for set_ in dataframe_to_update_[column_name_]:
                list_arr_.append(list(set_))
            dataframe_to_update_.update(pd.Series(list_arr_, name=column_name_))
'''


# In[ ]:





# ## DeepHGCal/TFRecords Data Preparation

# In[33]:


import tensorflow as tf


# In[34]:


MAX_RECHIT_LEN = 3600
MAX_TRACK_LEN = 120


# In[52]:


get_ipython().run_cell_magic(u'time', u'', u'import numpy as np\nfrom sklearn.preprocessing import StandardScaler\n\n\'\'\'This is the function to create graphs in the form readable by DeepHGCal\'\'\'\n\ndata_dict_list_ = []\nscaler = StandardScaler()\n\n# Global Features are track-based so they vary in length per-event\n# We find the maximum number of tracks that correspond to max_len of global feature vector\n# Is it a good idea to zero-pad global feature vectors less than max_len?\n\nfor event_id_ in range(number_of_events_):\n    data_dict_ = {}\n    senders_ = []\n    receivers_ = []\n    \n    track_event_df_ = track_global_df_[track_global_df_[\'event_id\'] == event_id_]\n    track_param_df_ = track_param_global_df_.loc[track_event_df_[\'track_id\']]\n    track_df_ = track_event_df_.merge(track_param_df_)\n\n    # Sort the tracks according to increasing track_eta and associate a label with each track\n    # This is done by resetting the track index based on increasing track_et\n    track_df_.sort_values(\'track_eta\', ascending=True, inplace=True)\n    track_df_.index = pd.RangeIndex(len(track_df_.index))  \n\n    rechit_event_df_ = rechit_global_df_[rechit_global_df_[\'event_id\']==event_id_]\n    rechit_param_df_ = rechit_param_global_df_[rechit_param_global_df_[\'event_id\']==event_id_]\n    rechit_param_df_.index = pd.RangeIndex(len(rechit_param_df_.index))  \n    if len(rechit_event_df_) != len(rechit_param_df_):\n        print("Error - param data and event data are not of equal length!")\n    \n    number_of_rechits_in_event_ = len(rechit_event_df_)\n    \n    # Set the node features as the track features that they belong to\n    node_indices_ = np.array(rechit_param_df_[\'rechit_local_id\'].tolist()).astype(int)\n    node_labels_ = []\n    \n    # Originally, we were setting node-level features based on the nodes but that can be done\n    # for the test data set; instead here we can set the node-level features as the track features\n    # at least for training and "learn" the track-level features (eta) based on which we can cluster the nodes?\n    # The question still remains how do we initialize the edges ???\n    \n    # Update: Reverting to node-level features for each node as of now \n    # Modify it to combine some form of track-level features (target?)\n    \n    # Skip the event if it has no rechits\n    if len(node_indices_) == 0:\n        print("Event ", event_id_, " has no rechits" )\n        continue\n        \n    rechit_feature_vector_ = np.transpose(np.array([\n        rechit_param_df_[\'rechit_r\'].tolist(),\n        rechit_param_df_[\'rechit_eta\'].tolist(),\n        rechit_param_df_[\'rechit_phi\'].tolist(),\n        rechit_param_df_[\'rechit_x\'].tolist(),\n        rechit_param_df_[\'rechit_y\'].tolist(),\n        rechit_param_df_[\'rechit_z\'].tolist(),]))\n    rechit_feature_vector_ = scaler.fit_transform(rechit_feature_vector_)\n    node_feature_vector_ = rechit_feature_vector_\n\n    # Initialize an array of zeros as labels for the nodes\n    # Replace each zero with the int (label) associated with the track to which the rechit belongs\n    # This becomes the target label to identify for that training example\n    node_label_array_ = [-1] * len(rechit_feature_vector_)\n    \n    # Associate a label with each rechit\n    for row in track_df_.itertuples():\n        # The itertuples() method for dataframes requires acccess by rows\n        # The matched_rechit_local_ids array is present in the sixth column of this \'row\'\n        # Do not confuse this with matched_rechit_global_ids that will cause a list index outofbounds error\n        track_rechit_id_array_ = row[6]\n                \n        # The index is the first element of this \'row\' and we use it as a unique label for the track\n        # Lower indices are thus associated with lower track eta\n        for id_ in track_rechit_id_array_:\n            node_label_array_[int(id_)] = row[0]\n\n    # Concatenate data and labels into a single matrix \n    # Note that we duplicate the labels to work with off-the-shelf DeepHGCal Model\n    nla_ = np.array(node_label_array_)\n    node_feature_vector_ = np.concatenate((rechit_feature_vector_, nla_[:, None], nla_[:, None]), axis = 1)\n\n    # data comprises of f features--here f is 6--and 1-d vector for labels (total f+1 columns)\n    # it has n rows of rechits making its shape n x f\n    # f is constant for all events but n varies with the event thus padding may be necessary\n    data_dict_ = {\n    "data": node_feature_vector_,\n    }\n    \n    data_dict_list_.append(data_dict_)\nprint(len(data_dict_list_), "graphs generated from data")')


# In[53]:


'''Convert the dict into a tf.Example then use a proto_buffer to 
serialize it into a compatible format for TFRecords'''

def create_tf_example(graph_dict=None, max_hits=None, data_dimensions=None):
    """
    :param graph_dict: dictionary with each key representing a feature vector
    :param labels: list with 
    Creates a tf.Example message ready to be written to a file.
    """
    if graph_dict is None:
        raise ValueError("No Rechit Feature Matrix provided")
    
    # Create a dictionary mapping the feature name to the tf.Example-compatible
    # data type.
    feature_matrix = {}
    # There is only one key in graph_dicts called 'data' but we keep this generalized
    key_count = 0
    for key in graph_dict.keys():
        # Add checks to see if additional keys are present (other than 'data')
        if key != 'data':
            if key_count == 0:
                raise ValueError("Key 'data' not found in rechit data dictionary")
            raise ValueError("Additional keys named ", key, " are present in data")
        # Presently, all features are concatenated into a matrix and converted to a tensor
        # We flatten this tensor and convert it into a FloatList that is then serialized
        # Pad the tensor with zeros to the MAX_HITS Margin
        if max_hits is None:
            MAX_HITS = 3600
        else:
            MAX_HITS = max_hits
        if data_dimensions is None:
            DATA_DIM = 6
        else:
            DATA_DIM = data_dimensions
        event_data_size = MAX_HITS * (DATA_DIM + 2)  # Add 2 columns to data dimensions to account for the label
        original_data = graph_dict[key].reshape(-1).astype(np.float32)
        
        # Pad the data with zeros to ensure the length of data from each event is constant
        # TODO: Can we improve this? After all, graphs are meant to be dynamic structures, right?
        if len(original_data) > event_data_size:
            padded_data = original_data[0:event_data_size]
        else:
            padded_data = np.concatenate((original_data, np.zeros(event_data_size - len(original_data))), axis=0)
        # Define the tf Feature to wrap the FloatList
        feature_matrix['data'] = tf.train.Feature(float_list=tf.train.FloatList(value=padded_data))
        
        # Increment key count to keep track of whether additional data is present 
        # This is used to throw a more informative error if key is not named 'data'
        key_count += 1
        
    # Create a Features message using tf.train.Example.
    example_proto = tf.train.Example(features=tf.train.Features(feature=feature_matrix))
    return example_proto.SerializeToString()


def _parse_function(example_proto):
    # Create a description of the features to be read from the TFRecord file(s).  
    feature_description = {
    'data': tf.FixedLenFeature([3200, 8], tf.float32),
    }

    # Parse the input tf.Example proto using the Feature dictionary above.
    return tf.parse_single_example(example_proto, feature_description)



# In[54]:


'''Write the TFRecord File'''

with tf.python_io.TFRecordWriter('tfrecords/ttbar-10.tfrecord', 
                                 options=tf.python_io.TFRecordOptions(
                                    tf.python_io.TFRecordCompressionType.GZIP)) as tfwriter:
    for event_number_, data_record_ in enumerate(data_dict_list_):
        # Undocumented function declaration; needs to be properly defined later
        tf_example_ = create_tf_example(data_record_, MAX_RECHIT_LEN, 8)
        tfwriter.write(tf_example_)


# In[ ]:




