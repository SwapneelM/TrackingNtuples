#!/usr/bin/env python
# coding: utf-8

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import uproot
import pandas as pd
from collections import OrderedDict

import argparse

parser = argparse.ArgumentParser(description='Parse the jobid for use in naming the outfiles', add_help=True)
parser.add_argument('jobid', type=int, help='The unique jobid that htcondor associates with each job')
parser.add_argument('n_events', type=int, help='The number of events contained in each root file')
parser.add_argument('outpath', type=str, help='The path to the output where you want to store your tfrecords')
args = parser.parse_args()

# ## View the Keys in the Imported Data

# In[111]:

jobid_ = args.jobid
number_of_events_ = args.n_events
outfile_ = jobid_ + ".root"
data_ = uproot.open(outfile_)["ntuples"]["tree"]
# data_.keys()


# ## Check the Integrity of the Imported Data 

# In[112]:


stereo_tp_idx_ = data_.array('stereoTPIndex')
mono_tp_idx_ = data_.array('monoTPIndex')
track_tp_idx_ = data_.array('trackTPIndex')

# Check that both have been generated for the same number of events
# Just for clarity
assert len(track_tp_idx_) == len(stereo_tp_idx_), "Track and Stereo Number of Events do not match"
assert len(track_tp_idx_) == len(mono_tp_idx_), "Track and Mono Number of Events do not match"
print ("\nTotal", len(track_tp_idx_), "events")


# # Optimisation Tests

# In[113]:


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


# In[114]:


mono_tp_idx_set_ = list_to_set(mono_tp_idx_)


# ## Load Data into Arrays

# In[115]:


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

# In[116]:


# Convert all tracking particle index lists to sets for faster search

mono_tp_idx_set_ = list_to_set(mono_tp_idx_)
stereo_tp_idx_set_ = list_to_set(stereo_tp_idx_)
track_tp_idx_set_ = list_to_set(track_tp_idx_)


# ## Preprocessing 2: Add all data into dataframes

# In[117]:


from pandas import DataFrame as df


# ### Create a Global Dataframe of Rechits

# In[118]:


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


# In[119]:


'''
Create the Global Rechit Array and Global Rechit Parameters Array'''
rechit_global_df_uncut_, rechit_param_global_df_uncut_ = create_global_rechit_df(
    stereo_tp_idx_, mono_tp_idx_, rechit_cartesian_, rechit_polar_)
#print rechit_global_df_.head(10)


# ## Place the Cuts (create DF for Graph Networks)

# In[120]:


'''Check the maximum number of hits in an event'''
max_len_ = 0 
for i in range(number_of_events_):
    len_idx_ = len(stereo_tp_idx_[i]) + len(mono_tp_idx_[i])
    if len_idx_ > max_len_:
        max_len_ = len_idx_
        #print max_len_
        
print ("Maximum hits in an event are: ", max_len_)
print("Average hits in an event are", float(sum([len(x) for x in stereo_tp_idx_]))/float(number_of_events_))     


# ### Format and Cut the Rechit DataFrames - also reorder Rechit Global and Local IDs

# In[121]:


'''This is done here to generate a reduced number of local indices for the tracks to match to rechits.
We will replace the rechit_global_df_ generated above and used below for matches with this new dataframe.'''

total_number_of_rechits_ = len(rechit_global_df_uncut_)
# Place the cuts on rechits not matched to any tracking particles (reduces rechits by about 75% for 10 events)
tp_unmatched_indices_ = rechit_global_df_uncut_[rechit_global_df_uncut_.rechit_tp_index.map(len) == 0].index
rechit_global_df_uncut_.drop(tp_unmatched_indices_, inplace=True)
rechit_param_global_df_uncut_.drop(tp_unmatched_indices_, inplace = True)
print("No of rechits left after 1st cut: ", len(rechit_global_df_uncut_))


# After each cut, we reset the global rechit ids so that 'iloc' and 'loc' do not throw errors
# Reset the Index of the Cut Dataframe that will become the new Global DataFrame
rechit_global_df_uncut_.index = pd.RangeIndex(len(rechit_global_df_uncut_.index))  
rechit_param_global_df_uncut_.index = pd.RangeIndex(len(rechit_global_df_uncut_.index))

# Update the Global Rechit IDs
rechit_global_id_dict_ = {}
rechit_global_id_dict_['rechit_id'] = range(len(rechit_global_df_uncut_))
rechit_global_df_uncut_.update(pd.DataFrame.from_dict(rechit_local_id_dict_))    
rechit_param_global_df_uncut_.update(pd.DataFrame.from_dict(rechit_local_id_dict_))    

# Place the cut on rechit eta to ensure you consider only the hits in the tracker
rechit_param_global_df_ = rechit_param_global_df_uncut_[np.abs(rechit_param_global_df_uncut_['rechit_eta']) <= 0.9].copy()
rechit_global_df_ = rechit_global_df_uncut_.iloc[rechit_param_global_df_['rechit_id'].index].copy()
print("No. of rechits left after second cut:", len(rechit_param_global_df_))

# After each cut, we reset the global rechit ids so that 'iloc' and 'loc' do not throw errors

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

print (len(rechit_param_global_df_), "of", total_number_of_rechits_, float(len(rechit_param_global_df_))/float(total_number_of_rechits_), "hits remain")


# ## Match the Rechits to Tracks and Create a Global Array of Tracks

# In[122]:


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

# In[123]:


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

# In[124]:


from mpl_toolkits.mplot3d import Axes3D
from cycler import cycler
from matplotlib.colors import Colormap

#fig_ = plt.figure()
#ax_ = Axes3D(fig_)


# ### Analyse Matched/Unmatched Rechits

# In[125]:


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

# In[126]:


track_count_, track_ids_ = count_matched_items('track')
rechit_count_, rechit_ids_ = count_matched_items('rechit')


# ## DeepHGCal/TFRecords Data Preparation

# In[130]:


import tensorflow as tf


# In[131]:


max_hits = 3600
max_tracks = 100


# In[132]:


from sklearn.preprocessing import StandardScaler, OneHotEncoder

'''This is the function to create graphs in the form readable by DeepHGCal'''

data_dict_list_ = []
scaler = StandardScaler()

# Global Features are track-based so they vary in length per-event
# We find the maximum number of tracks that correspond to max_len of global feature vector
# Is it a good idea to zero-pad global feature vectors less than max_len?

for event_id_ in range(number_of_events_):
    data_dict_ = {}
    senders_ = []
    receivers_ = []
    
    track_event_df_ = track_global_df_[track_global_df_['event_id'] == event_id_]
    track_param_df_ = track_param_global_df_.loc[track_event_df_['track_id']]
    track_df_ = track_event_df_.merge(track_param_df_)

    # Sort the tracks according to increasing track_eta and associate a label with each track
    # This is done by resetting the track index based on increasing track_et
    track_df_.sort_values('track_eta', ascending=True, inplace=True)
    track_df_.index = pd.RangeIndex(len(track_df_.index))  

    rechit_event_df_ = rechit_global_df_[rechit_global_df_['event_id']==event_id_]
    rechit_param_df_ = rechit_param_global_df_[rechit_param_global_df_['event_id']==event_id_]
    rechit_param_df_.index = pd.RangeIndex(len(rechit_param_df_.index))  
    if len(rechit_event_df_) != len(rechit_param_df_):
        print("Error - param data and event data are not of equal length!")
    
    number_of_rechits_in_event_ = len(rechit_event_df_)
    
    # Set the node features as the track features that they belong to
    node_indices_ = np.array(rechit_param_df_['rechit_local_id'].tolist()).astype(int)
    node_labels_ = []
    
    # Originally, we were setting node-level features based on the nodes but that can be done
    # for the test data set; instead here we can set the node-level features as the track features
    # at least for training and "learn" the track-level features (eta) based on which we can cluster the nodes?
    # The question still remains how do we initialize the edges ???
    
    # Update: Reverting to node-level features for each node as of now 
    # Modify it to combine some form of track-level features (target?)
    
    # Skip the event if it has no rechits
    if len(node_indices_) == 0:
        print("Event ", event_id_, " has no rechits" )
        continue
        
    rechit_feature_vector_ = np.transpose(np.array([
        rechit_param_df_['rechit_r'].tolist(),
        rechit_param_df_['rechit_eta'].tolist(),
        rechit_param_df_['rechit_phi'].tolist(),
        rechit_param_df_['rechit_x'].tolist(),
        rechit_param_df_['rechit_y'].tolist(),
        rechit_param_df_['rechit_z'].tolist(),
        np.zeros(len(rechit_param_df_['rechit_eta'])),
        np.zeros(len(rechit_param_df_['rechit_eta'])),
        np.zeros(len(rechit_param_df_['rechit_eta'])),
        np.zeros(len(rechit_param_df_['rechit_eta'])),
    ]))
    
    track_feature_vector_ = np.transpose(np.array([
        np.zeros(len(track_df_['track_eta'])),
        track_df_['track_eta'].tolist(),
        track_df_['track_phi'].tolist(),
        np.zeros(len(track_df_['track_eta'])),
        np.zeros(len(track_df_['track_eta'])),
        np.zeros(len(track_df_['track_eta'])),
        track_df_['track_dsz'].tolist(),
        track_df_['track_dxy'].tolist(),
        track_df_['track_qoverp'].tolist(),
        track_df_['track_pt'].tolist(),
    ]))

    node_feature_vector_ = np.vstack((rechit_feature_vector_, track_feature_vector_))
    node_feature_vector_ = scaler.fit_transform(node_feature_vector_)    
    
    # Initialize an array of zeros as labels for the nodes
    # Replace each zero with the int (label) associated with the track to which the rechit belongs
    # This becomes the target label to identify for that training example
    # Later, we will one_hot_encode the target label for working with the tensorflow model
    node_label_array_ = [0] * len(rechit_feature_vector_)
    track_labels_ = []
    
    # Associate a label with each rechit
    for trk_idx_, row in track_df_.iterrows():
        track_rechit_id_array_ = row.rechit_local_ids
        # We use the rechit local index as a unique label for the track
        # Lower indices are thus associated with lower track eta
        for id_ in track_rechit_id_array_:
            node_label_array_[int(id_)] = trk_idx_ + 1  # We can now use index '0' as a 'noise' label
        
        # Also associate a label with each track that we consider as a point among the data
        track_labels_.append(trk_idx_ + 1)
    
    print(len(node_label_array_))   
    #print(node_indices_.shape)
    #node_labels_ = np.vstack((node_indices_, label_array_)).T
    #print(node_labels_.shape)
    #assert (len(node_labels_)==len(rechit_feature_vector_)), "Node label and Node feature vector length mismatch"
    
    # Concatenate node and track labels into a single array
    # Thereafter, concatenate everything into the node feature matrix 
    # Note that we duplicate the labels to work with off-the-shelf DeepHGCal Model
    nla_ = np.concatenate(((np.array(node_label_array_)), np.array(track_labels_)), axis=0)  # add track labels here
    # node_feature_vector_ = np.concatenate((rechit_feature_vector_, nla_[:, None], nla_[:, None]), axis = 1)
    assert node_feature_vector_.shape[0] == nla_.shape[0], "Labels are not equal to training nodes/examples"
    # data comprises of f features--here f is 6--and 1-d vector for labels (total f+1 columns)
    # it has n rows of rechits making its shape n x f
    # f is constant for all events but n varies with the event thus padding may be necessary
    data_dict_ = {
    "data": node_feature_vector_,
    "labels": nla_,
    }
    
    data_dict_list_.append(data_dict_)
print(len(data_dict_list_), "graphs generated from data")


# In[135]:


from sklearn.preprocessing import OneHotEncoder

'''Convert the dict into a tf.Example then use a proto_buffer to 
serialize it into a compatible format for TFRecords'''

def create_tf_example(graph_dict=None, max_hits=None, max_tracks=None, set_one_hot_labels=True):
    """
    :param graph_dict: dictionary with each key representing a feature vector
    :param labels: list with 
    Creates a tf.Example message ready to be written to a file.
    """
    one_hot_encoder = OneHotEncoder(handle_unknown='ignore')
    
    if max_hits is None:
        max_hits = 3600
    
    if graph_dict is None:
        raise ValueError("No Rechit Feature Matrix provided")

    # TODO: Can we eliminate this padding? After all, graphs are meant to be dynamic structures, right?
    # Presently, all features are concatenated into a matrix and converted to a tensor
    if 'data' in graph_dict.keys():
        original_data = graph_dict['data']
        data_dimensions = original_data.shape[1]
        # If number of points are less than max_hits threshold
        if original_data.shape[0] < max_hits:
            event_data_size = max_hits * data_dimensions
            
            # Reshape the original data into a 1D array
            # TODO: Is np.ravel faster? Does it matter?
            original_data = graph_dict['data'].reshape(-1).astype(np.float32)
            # Pad the data with zeros to ensure the length of data from each event is constant
            padded_data = np.concatenate((original_data, np.zeros(event_data_size - len(original_data))), axis=0)
            padded_data = np.reshape(padded_data, (max_hits, data_dimensions))
        else:
            # Eliminate the first bits of the data (these are definitely rechits)
            # The last part of the data is the tracks and we do not want to eliminate those
            # Since they are useful as 'centroids' for the clustering
            padded_data = original_data[-event_data_size:]
    else:
        raise ValueError("Key 'data' not found in rechit data dictionary")
    
    # Handle the labels either separately as one-hot-encoded data or together
    if 'labels' in graph_dict.keys():
        node_labels = graph_dict['labels']
        
        # Zero-pad the labels (in the case of one-hot encoded labels, the '0' label
        # corresponds to noisy points in the cloud)
        if len(node_labels) > max_hits:
            padded_labels = node_labels[-max_hits:]
        else:
            # Pad the data with zeros to ensure the length of data from each event is constant
            padded_labels = np.concatenate((node_labels, np.zeros(max_hits - len(node_labels))))
        
        # Check if one-hot-encoded labels are desired
        if set_one_hot_labels is True:
            if max_tracks is None:
                # Set default value of maximum tracks
                # This will add as many dimensions to your data
                max_tracks = 100
            
            # Create one-hot-encoded representation of the target labels for the data
            target_labels = one_hot_encoder.fit_transform(np.array([padded_labels]).T).toarray()
            # If there are more labels than max_tracks, just cut out the last set of labels
            if target_labels.shape[1] > max_tracks:
                target_labels = target_labels[:, :max_tracks]
            else:
                # Append zero-columns as padding to these data items
                padding_columns = np.zeros((max_hits, max_tracks-target_labels.shape[1]))
                target_labels = np.hstack((target_labels, padding_columns))
        else:            
            # Code to store regular labels in a duplicated array for DeepHGCal
            target_labels = np.concatenate((padded_labels[:, None], padded_labels[:, None]), axis = 1)
    else:
        raise ValueError("Key 'labels' not found in rechit data dictionary")
    
    # Once the labels and data is padded and concatenated, finalize the data into a 2D matrix
    final_data = np.hstack((padded_data, target_labels))
    print(final_data.shape)
    # Create a dictionary mapping the feature name to the tf.Example-compatible
    # data type.
    feature_matrix = {}
        
    # We flatten this tensor and convert it into a FloatList that is then serialized
    # Define the tf Feature to wrap the FloatList
    feature_matrix['data'] = tf.train.Feature(float_list=tf.train.FloatList(value=final_data.ravel()))
        
    # Create a Features message using tf.train.Example.
    example_proto = tf.train.Example(features=tf.train.Features(feature=feature_matrix))
    return example_proto.SerializeToString()


def _parse_function(example_proto, data_dimensions=None):
    '''If you would like to read and parse the data stored in TFRecord format, 
    refer to the [Dev] Prototyping Graph Neural Networks Notebook'''
    # max_tracks = 100; features = 10
    if data_dimensions is None:
        data_dimensions = (3600, 110)
        
    # Create a description of the features to be read from the TFRecord file(s).  
    feature_description = {
    'data': tf.FixedLenFeature(data_dimensions, tf.float32),
    }

    # Parse the input tf.Example proto using the Feature dictionary above.
    return tf.parse_single_example(example_proto, feature_description)



# In[136]:
import os

output_directory = os.path.join(args.outpath, 'tfrecords')

if not os.path.exists(output_directory):
    os.makedirs(output_directory)

'''Write the TFRecord File'''

with tf.python_io.TFRecordWriter(output_directory + '/' + str(args.jobid) + '_' + str(args.n_events) + '.tfrecord', 
                                 options=tf.python_io.TFRecordOptions(
                                    tf.python_io.TFRecordCompressionType.GZIP)) as tfwriter:
    for event_number_, data_record_ in enumerate(data_dict_list_):
        tf_example_ = create_tf_example(data_record_)
        tfwriter.write(tf_example_)


# In[ ]:




