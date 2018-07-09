# MEG contrast discrimination script
# Last edited 25/10/16

# INFO:
# The script is set up so you can quit at any time by pressing ESCAPE and then waiting for the next trial
# HOWEVER, when restarting the experiment, will need to change the number of blocks in the GUI (defaults is 4)
# Presents four identical blocks of randomised (beforehand) trials
# ISI and ITI are randomised
# A condition file is saved first so that if the experiment is quit and restarted, the same conditions can be used in the next blocks
# A data file for each block is saved
# And also a trigger file is saved for the whole session
# If the experiment is restarted, the new trigger file is saved with 'newsession' in the name of the file

# to do:
# - parallel port stuff

from __future__ import division  # so that 1/3=0.333 instead of 1/3=0
from psychopy import locale_setup, visual, core, data, event, logging, gui, monitors

# sound,
import parallel #NB not psychopy's parallel module -- FOR TRIGGERS AND RESPONSES
from psychopy.constants import *  # things like STARTED, FINISHED
import numpy as np  # whole numpy lib is available, prepend 'np.'
from numpy import sin, cos, tan, log, log10, pi, average, sqrt, std, deg2rad, rad2deg, linspace, asarray
from numpy.random import random, randint, normal, shuffle
import os  # handy system and path functions
import sys # to get file system encoding
import random
import csv
import time as trackTime

#print monitors.calibTools.monitorFolder

## -- NB hard code monitor settings
monitors.calibTools.monitorFolder = '/groups/Projects/P1314/MEG/monitors/'
##
#print monitors.calibTools.monitorFolder


# GUI for inputting subject data:
subjectDetails = {}
subjectDetails['subjectNumber'] = ''
subjectDetails['gender'] = ''
subjectDetails['age'] = ''
subjectDetails['blocks'] ='4'

dlg = gui.DlgFromDict(subjectDetails, title='SubjectInfo', order=['subjectNumber','gender','age', 'blocks'])

if not dlg.OK:
    print('User cancelled the experiment')
    core.quit()

# SETTINGS *****************************
dopracticetrials = True
blocks = int(subjectDetails['blocks'])
reps = 160 # 320 of each contrast condition in total (160 per interval configuration)
totalBlockTrials = int(reps*6/blocks) # 240 (with grand total of 960)

MonFrameRate = 60.0 # give monitor refresh rate
MonSizeHor = 53.0 # centimeters horizontally
MonSizeVer = 41.0 # centimeters vertically
MonPixels = [1024.0, 768.0]
MonDist = 105.0 # participant distance from monitor

# stimulus specific settings
gratingSF = 10.0 # cycles in the stimulus (right now 10 cycles in a 10deg stimulus)
gratingSize = 10.0 # visual angle degrees

# calculating pixels per visual angle, ect
tempRes1 = (MonPixels[0])
tempRes2 = (MonPixels[1])

print tempRes1
print MonSizeHor
print MonDist

PixVisAng =  35#for MEG, at viewdist = 105cm
#PixVisAng2 = 35#for MEG, at viewdist = 105cm

# checking the calculated pixels per visual angle are consistent horizontally and vertically
#if PixVisAng != PixVisAng2:
#    print "%%%%%%%%%%%% WARNING!!! %%%%%%%%%%%%"
#    print "     Pixel size is inconsistent     "
#    print "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"
#    core.quit()

# Checking if we already have a condition file (this will happen if we unexpectedly quit the experiment previously)
conditionfile = './Results/S%s_conditions.npy' % (subjectDetails['subjectNumber'])

if os.path.isfile(conditionfile):
    a = np.load(conditionfile)
    contrast_g1 = a[0]
    contrast_g2 = a[1]
    correctResps = a[2]
    indexing = a[3]
else:
    # creating conditions
    contrast_g1 = [0.5, 0.54, 0.66, 0.5, 0.5, 0.5]
    contrast_g2 = [0.5, 0.5, 0.5, 0.5, 0.54, 0.66]
    
    correctResps = ["1", "1", "1", "2", "2", "2"]
    
    condCodes_g1 = [11, 21, 31, 11, 11, 11]
    condCodes_g2 = [11, 11, 11, 11, 21, 31]
    
    nooftiles = int(totalBlockTrials/len(correctResps))
    
    contrast_g1 = np.tile(contrast_g1,nooftiles)
    contrast_g2 = np.tile(contrast_g2,nooftiles)
    
    correctResps = np.tile(correctResps,nooftiles)
    
    condCodes_g1 = np.tile(condCodes_g1,nooftiles)
    condCodes_g2 = np.tile(condCodes_g2,nooftiles)
    
    indexing = np.random.permutation(range(int(totalBlockTrials)))
    
    # saving the conditions for later use in case of problem during experiment:
    conditionInfo = [] # making a list of four arrays
    conditionInfo.append(contrast_g1)
    conditionInfo.append(contrast_g2)
    conditionInfo.append(correctResps)
    conditionInfo.append(indexing)
    np.save(conditionfile, conditionInfo)

# making a window
win = visual.Window(
    size=MonPixels,
    units="pix",
    fullscr=True,
    allowGUI=False,
    monitor='MEG_Screen'
)

print (PixVisAng)
print (gratingSize)
print (gratingSF)

print gratingSize/(PixVisAng*gratingSF)

# making the gratings
grating1 = visual.GratingStim(
    win=win,
    size=[gratingSize*PixVisAng, gratingSize*PixVisAng],
    mask="raisedCos",
    maskParams = {'fringeWidth':0.2},
    units="pixels",
    ori=90,
    sf = gratingSize/(PixVisAng*gratingSF) # 1/((PixVisAng*gratingSF)/gratingSize)
)

grating2 = visual.GratingStim(
    win=win,
    size=[gratingSize*PixVisAng, gratingSize*PixVisAng],
    mask="raisedCos",
    maskParams= {'fringeWidth':0.2},
    units="pixels",
    ori=90,
    sf = gratingSize/(PixVisAng*gratingSF)
)


# Stimulus and interval length:
ISI = int(0.5*MonFrameRate)
ITI = 1.1
stimLength = round(0.1 * MonFrameRate) # 100ms times refresh rate equals the number of frames to present stimulus for
stimLength = int(stimLength)


# make/show the welcome screen + practice trials here
posint = int(tempRes2/8) # an eighth of the number of vertical pixels; used for positioning text on screen

welcometext = visual.TextStim(
win=win,
color=[1, 1, 1],
pos = [0, posint*2.5],
)


### --- set up trigger port
# AG code -- function to send a trigger
parport = parallel.Parallel('/dev/parport0')  # address for parallel port on many machines
#MUST set the read/write mode in linux, 0=read 1=write 
parport.setDataDir(1)
#set the parallel port data pins (2-9) to zero before we start
parport.setData(0)


### ---- set up reponse port
responsePort = parallel.Parallel('/dev/parport1') # right hand lumitouch response box
#MUST set the read/write mode in linux, 0=read 1=write 
responsePort.setDataDir(0)


def getLumitouchResponseMEG():
    while int(responsePort.getData()) == 0:
        trackTime.sleep(0.001)
        
        quitExp = event.getKeys(keyList=['escape'])
        if quitExp:
            print('User has exited the experiment')
            win.close()
            core.quit()
    
    responseReceived = str(int(responsePort.getData())) # code currently expects a string
    return responseReceived
        

def sendTriggerMEG(triggerValue):
    #this should be called immediately after win.flip for any key stimulus
    # takes an integer argument for the trigger value to be sent
    parport.setData(int(triggerValue))
### --------------




# practice trials:
if dopracticetrials == True:
    
    welcometext.text='Welcome to our MEG experiment!'
    welcometext.draw()
    welcometext.text='During the experiment you will see two gratings appear one after the other in the centre of the screen. You will need to choose the grating that is HIGHER in contrast.'
    welcometext.pos= [0, posint]
    welcometext.draw()
    welcometext.text='Please wait until you have seen both gratings and then press LEFT button for the FIRST grating or RIGHT button for the SECOND grating.'
    welcometext.pos= [0, 0]
    welcometext.draw()
    welcometext.text='Before the experiment begins you will have 5 practice trials with feedback (correct/incorrect). Note that you will not have this feedback during the actual experiment.'
    welcometext.pos= [0, -(posint)]
    welcometext.draw()
    welcometext.text="We'll start some practice trials for you shortly (Enter)"
    welcometext.pos= [0, -(posint*2)]
    welcometext.draw()
    win.flip()
    event.waitKeys()

    welcometext.text="Starting practice trials .... (Enter)"
    welcometext.pos= [0,0]
    welcometext.draw()
    win.flip()
    core.wait(3.0)

    win.flip()
    core.wait(1.5)

    practiceContrast1 = [0.72, 0.5, 0.64, 0.58, 0.5]
    practiceContrast2 = [0.5, 0.68, 0.5, 0.5, 0.66]
    practiceCorr = ['1', '2', '1', '1', '2']
    
    
    for practice in range(5):
        # setting up grating contrast# to do trigger levels:
        grating1.contrast = practiceContrast1[practice]
        grating2.contrast = practiceContrast2[practice]
        correctResp = practiceCorr[practice]
        # showing stimuli:
        for run1 in range(stimLength):
            grating1.draw()
            win.flip()
        
        for run2 in range(ISI):
            win.flip()
            
        for run3 in range(stimLength):
            grating2.draw()
            win.flip() # makes second grating appear
        win.flip() 
        clock = core.Clock() # make a clock
        ##practkeypress = event.waitKeys(keyList=['1','2'],timeStamped=clock) # waiting for key to be pressed
        ##practkeytime = practkeypress[0][1]
        ##practkeytime = round(practkeytime, 3)
        ##practrespKey = practkeypress[0][0]
        
        practrespKey = getLumitouchResponseMEG()
        practkeytime = clock.getTime()
        
        # make sure it waits the correct amount of time before starting a new trial:
        if practkeytime < ITI:
            timewait = ITI - practkeytime
            core.wait(timewait)
        
        # correct or incorrect?
        if practrespKey == correctResp:
            practtext = 'Correct!'
            practcol = [-1, 0, -1]
        else:
            practtext = 'Incorrect'
            practcol = [0, -1, -1]
        
        win.flip()
        welcometext.text = practtext
        welcometext.color = practcol
        welcometext.pos = [0, 0]
        welcometext.draw()
        win.flip()
        core.wait(1)
    win.flip()
    core.wait(1)

# end of practice screen
welcometext.text = 'End of practice. The experiment is wil start soon. Remember:'
welcometext.color = [1, 1, 1]
welcometext.pos = [0, posint]
welcometext.draw()
welcometext.text = '- There will be no feedback and so the trials will come in quick succession'
welcometext.pos = [0, 0]
welcometext.draw()
welcometext.text = '- Press LEFT for FIRST grating'
welcometext.pos = [0, -(posint/2)]
welcometext.draw()
welcometext.text = '- Press RIGHT for SECOND grating'
welcometext.pos = [0, -(posint)]
welcometext.draw()
welcometext.text = '**** We will start the experiment for you soon ... (Enter)****'
welcometext.pos = [0, -(posint*3)]
welcometext.draw()
win.flip()
event.waitKeys()
win.flip()
core.wait(1.5)
win.flip()

welcometext.pos = [0,0]
welcometext.text = '**** Get ready to start ****'
welcometext.draw()
win.flip()

core.wait(2.5)
win.flip()
core.wait(2.5)

if (blocks == 4):
    triggerfile = 'Results/S%s_MEGtriggers.csv' % (subjectDetails['subjectNumber'])
elif (blocks < 4):
    triggerfile = 'Results/S%s_MEGtriggers_newsession.csv' % (subjectDetails['subjectNumber'])
# making the trigger file:
t = open(triggerfile, 'w')
# Writing column headings to the trigger file
t.write('Subject number, %s\n\n\n' % (subjectDetails['subjectNumber']))
t.write('Trigger code, Timing\n')



# start of block loop
for bl in range(blocks):
    print bl
    if (blocks == 4):
        # making a normal new data file for saving conditions and responses; also a separate file for triggers
        datafile = 'Results/S%s_contrastMEG_%s.csv' % (subjectDetails['subjectNumber'], str(bl))
    elif (blocks < 4):
        # correcting block numbers if we are running less than 4 blocks:
        blockcorrect = 4-blocks
        datafile = 'Results/S%s_contrastMEG_%s.csv' % (subjectDetails['subjectNumber'], str((bl+blockcorrect)))
    else:
        print 'Max number of blocks exceded'
        core.quit()
    
    
    
    print datafile
    
    
    # Writing column headings to the data file
    d = open(datafile, 'w')
    
    d.write('Subject number, %s\n' % (subjectDetails['subjectNumber']))
    d.write('Gender, %s\n' % (subjectDetails['gender']))
    d.write('Age, %s\n\n\n' % (subjectDetails['age']))
    d.write('Trial,Grating1,Grating2,CorrResp,Resp,IsCorr,RT\n')
    
    #-----------------send start of block trigger here-----------------
    expclock = core.Clock() # make a clock
    trigger = 1
    sendTriggerMEG(trigger)
    time = expclock.getTime()
    t.write('%d, %3.6f\n' % (trigger, time))
    # reset trigger after 
    win.flip()
    trigger = 0
    sendTriggerMEG(trigger)
    #------------------------------------------------------------------
    
    # trial loop
    for trial in range(totalBlockTrials): # range(totalBlockTrials*blocks) # minus the placeholder so that we run fewer trials when we have already done some
        # length of presentations:
        ISI = round((0.6 + random.randrange(0, 1000)/5000) * MonFrameRate) # picks a random length (0.6 - 0.8s) with millisecond precision
        ITI = round((1 + random.randrange(0, 1000)/5000)) # picks a random length (1 - 1.2s) with millisecond precision
        # ITI is in seconds and not frames because easier to wait for keys
        
        ISI = int(ISI)
        
        # setting up grating contrast levels:
        g1 = contrast_g1[indexing[trial]]
        g2 = contrast_g2[indexing[trial]]
        grating1.contrast = g1
        grating2.contrast = g2
        correctResp = correctResps[indexing[trial]]
        # showing stimuli:
        #-----------------send first interval trigger here-----------------
        trigger = condCodes_g1[indexing[trial]]
        #------------------------------------------------------------------
        curCount=0
        for run1 in range(stimLength):
            grating1.draw()
            win.flip()
            if curCount == 0:
                time = expclock.getTime()
                t.write('%d, %3.6f\n' % (trigger, time))
                sendTriggerMEG(trigger) #trigger code
                curCount=1
                
        trigger = 0
        curCount=0
        for run2 in range(ISI):
            win.flip()
            if curCount == 0:
                sendTriggerMEG(trigger) #trigger code
                curCount=1
            
        #-----------------send second interval trigger here-----------------
        trigger = condCodes_g2[indexing[trial]]

        #------------------------------------------------------------------
        curCount=0
        for run3 in range(stimLength):
            grating2.draw()
            win.flip() # makes second grating appear
            if curCount == 0:
                time = expclock.getTime()
                t.write('%d, %3.6f\n' % (trigger, time))
                sendTriggerMEG(trigger) #trigger code
                curCount=1
        
        win.flip()
        
        #reset trigger to zero
        trigger = 0
        sendTriggerMEG(trigger) #trigger code
        
        clock = core.Clock() # make a clock
        ##keypress = event.waitKeys(keyList=['1','2'],timeStamped=clock) # waiting for key to be pressed
        ##keytime = keypress[0][1]
        ##keytime = round(keytime, 6)
        ##respKey = keypress[0][0]
        
        respKey = getLumitouchResponseMEG()
        keytime = clock.getTime()
        
        # correct or incorrect?
        if respKey == correctResp:
            correct = 1
            trigger = 101 # correct
        else:
            correct = 0
            trigger = 201 # incorrect
        
        ### log event as a trigger
        ##sendTriggerMEG(trigger) #trigger code
        
        time = expclock.getTime()
        t.write('%d, %3.6f\n' % (trigger, time))
        # make sure it waits the correct amount of time before starting a new trial:
        if keytime < ITI:
            timewait = ITI - keytime
            core.wait(timewait)
        
        
        # if we click escape, quit the experiment:
        quitExp = event.getKeys(keyList=['escape'])
        if quitExp:
            print('User has exited the experiment')
            win.close()
            core.quit()
        
        # saving the data into our csv file
        d.write('%d, %.2f, %.2f, %s, %s, %d, %1.3f\n' % ((trial+1), g1, g2, correctResp, respKey, correct, keytime))
        
        # if all the trials in the block are done, do this:
        if (trial == (totalBlockTrials-1)):
             # if all the blocks are done too, quit:
            if (bl == (blocks-1)):
                text = visual.TextStim(win, text='End of experiment\n\nPLEASE KEEP STILL UNTIL WE TELL YOU TO MOVE.\n\n Thank you. (Enter to end)')
                d.close()
                t.close()
                text.draw()
                win.flip()
                event.waitKeys()
                win.close()
            else:
                # if not, move on to next block:
                text = visual.TextStim(win, text='End of block. Click to start a new block once you have rested')
                d.close()
                text.draw()
                win.flip()
                ##event.waitKeys()
                respKey = getLumitouchResponseMEG()
                win.flip()


