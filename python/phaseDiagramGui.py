import PySimpleGUI as sg
import json
import matplotlib.pyplot as plt
import numpy as np
import math
import os
import subprocess

def processPhaseDiagramData(fname, elx, ts, x1, x2, p1, p2, mint, maxt):
    f = open(fname,)
    data = json.load(f)
    f.close()
    if list(data.keys())[0] != '1':
        print('Output does not contain data series')
        exit()
    for i in list(data.keys()):
        mint = min(mint,data[i]['temperature'])
        maxt = max(maxt,data[i]['temperature'])
        if (data[i]['# solution phases'] + data[i]['# pure condensed phases']) == 2:
            ts.append(data[i]['temperature'])
            boundPhases = []
            boundComps = []
            for phaseName in list(data[i]['solution phases'].keys()):
                if (data[i]['solution phases'][phaseName]['moles'] > 0):
                    boundPhases.append(phaseName)
                    boundComps.append(data[i]['solution phases'][phaseName]['elements'][elx]['mole fraction of phase by element'])
            for phaseName in list(data[i]['pure condensed phases'].keys()):
                if (data[i]['pure condensed phases'][phaseName]['moles'] > 0):
                    boundPhases.append(phaseName)
                    boundComps.append(data[i]['pure condensed phases'][phaseName]['elements'][elx]['mole fraction of phase by element'])
            x1.append(boundComps[0])
            x2.append(boundComps[1])
            p1.append(boundPhases[0])
            p2.append(boundPhases[1])
    return mint, maxt

def runCalc(ts, x1, x2, p1, p2, mint, maxt):
    print('Thermochimica calculation initiated.')
    subprocess.run(['./bin/PhaseDiagramDataGen',filename])
    print('Thermochimica calculation finished.')

    fname = 'thermoout.json'

    mint, maxt = processPhaseDiagramData(fname, el2, ts, x1, x2, p1, p2, mint, maxt)

    boundaries = []
    b = []
    for i in range(len(p1)):
        # If a miscibility gap label has been used unnecessarily, remove it
        if p1[i].find('#2') > 0:
            if not(p1[i][0:p1[i].find('#2')] == p2[i]):
                p1[i] = p1[i][0:p1[i].find('#2')]
        if p2[i].find('#2') > 0:
            if not(p2[i][0:p2[i].find('#2')] == p1[i]):
                p2[i] = p2[i][0:p2[i].find('#2')]

        repeat = False
        for j in range(len(boundaries)):
            if (boundaries[j][0] == p1[i]) and (boundaries[j][1] == p2[i]):
                b.append(j)
                repeat = True
        if not(repeat):
            boundaries.append([p1[i],p2[i]])
            b.append(len(boundaries)-1)

    # Start figure
    fig = plt.figure()
    ax = fig.add_axes([0.2, 0.1, 0.75, 0.85])
    for j in range(len(boundaries)):
        inds = [i for i, k in enumerate(b) if k == j]
        ax.plot(np.array(x1)[inds],np.array(ts)[inds],'.')
        ax.plot(np.array(x2)[inds],np.array(ts)[inds],'.')
        minj = np.argmin(np.array(ts)[inds])
        maxj = np.argmax(np.array(ts)[inds])
        if (np.array(ts)[inds][minj] > mint):
            ax.plot([np.array(x1)[inds][minj],np.array(x2)[inds][minj]],[np.array(ts)[inds][minj],np.array(ts)[inds][minj]],'k-')
        if (np.array(ts)[inds][maxj] < maxt):
            ax.plot([np.array(x1)[inds][maxj],np.array(x2)[inds][maxj]],[np.array(ts)[inds][maxj],np.array(ts)[inds][maxj]],'k-')
    ax.set_xlim(0,1)
    plt.show()
    return mint, maxt

def writeInputFile(filename,xlo,xhi,tlo,thi,pressure,tunit,punit,munit,el1,el2,datafile):
    with open(filename, 'w') as inputFile:
        inputFile.write('! Python-generated input file for Thermochimica\n')
        xstep = (float(xhi)-float(xlo))/float(nxstep)
        inputFile.write('x          = ' + str(xlo) + ':' + str(xhi) + ':' + str(xstep) + '\n')
        tstep = (float(thi)-float(tlo))/float(ntstep)
        inputFile.write('temperature          = ' + str(tlo) + ':' + str(thi) + ':' + str(tstep) + '\n')
        inputFile.write('pressure          = ' + str(pressure) + '\n')
        inputFile.write('temperature unit         = ' + tunit + '\n')
        inputFile.write('pressure unit          = ' + punit + '\n')
        inputFile.write('mass unit          = \'' + munit + '\'\n')
        inputFile.write('iEl         = ' + str(atomic_number_map.index(el1)+1) + ' ' + str(atomic_number_map.index(el2)+1) + '\n')
        inputFile.write('data file         = ' + datafile + '\n')

atomic_number_map = [
    'H','He','Li','Be','B','C','N','O','F','Ne','Na','Mg','Al','Si','P',
    'S','Cl','Ar','K','Ca','Sc','Ti','V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
    'Ga','Ge','As','Se','Br','Kr','Rb','Sr','Y','Zr','Nb','Mo','Tc','Ru','Rh',
    'Pd','Ag','Cd','In','Sn','Sb','Te','I','Xe','Cs','Ba','La','Ce','Pr','Nd',
    'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb','Lu','Hf','Ta','W','Re',
    'Os','Ir','Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
    'Pa','U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm','Md','No','Lr','Rf','Db',
    'Sg','Bh','Hs','Mt','Ds','Rg','Cn','Nh','Fl','Mc','Lv','Ts', 'Og'
]

file_list_column = [
    [
        sg.Text("Database Folder"),
        sg.In(size=(25, 1), enable_events=True, key="-FOLDER-"),
        sg.FolderBrowse(),
    ],
    [
        sg.Listbox(
            values=[], enable_events=True, size=(40, 20), key="-FILE LIST-"
        )
    ],
]

dataWindow = sg.Window('Thermochimica database selection', file_list_column, location = [0,0], finalize=True)
folder = os.getcwd()+'/data'
try:
    file_list = os.listdir(folder)
except:
    file_list = []

fnames = [
    f
    for f in file_list
    if os.path.isfile(os.path.join(folder, f))
    and f.lower().endswith((".dat", ".DAT"))
]
dataWindow["-FILE LIST-"].update(fnames)

timeout = 50
inputSize = 20

while True:
    event, values = dataWindow.read()
    if event == sg.WIN_CLOSED or event == 'Exit':
        break
    elif event == "-FOLDER-":
        folder = values["-FOLDER-"]
        try:
            file_list = os.listdir(folder)
        except:
            file_list = []
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(folder, f))
            and f.lower().endswith((".dat", ".DAT"))
        ]
        dataWindow["-FILE LIST-"].update(fnames)
    elif event == "-FILE LIST-":  # A file was chosen from the listbox
        try:
            datafile = os.path.join(
                folder, values["-FILE LIST-"][0]
            )
            with open(datafile) as f:
                f.readline() # read comment line
                line = f.readline() # read first data line (# elements, # phases, n*# species)
                nElements = int(line[1:5])
                nSoln = int(line[6:10])
                elements = []
                for i in range(math.ceil((nSoln+3)/15)-1):
                    f.readline() # read the rest of the # species but don't need them)
                for i in range(math.ceil(nElements/3)):
                    els = f.readline() # read a line of elements (3 per line)
                    elLen = 25 # formatted 25 wide
                    for j in range(3):
                        elements.append(els[1+j*elLen:(1+j)*elLen].strip())
            i = 0
            while i < nElements:
                try:
                    index = atomic_number_map.index(elements[i])+1 # get element indices in PT (i.e. # of protons)
                    i = i + 1
                except ValueError:
                    print(elements[i]+' not in list') # if the name is bogus (or e(phase)), discard
                    elements.remove(elements[i])
                    nElements = nElements - 1
            elSelectLayout = [sg.Column([[sg.Text('Element 1')],[sg.Combo(elements[:nElements],default_value=elements[0],key='-el1-')]],vertical_alignment='t'),
                              sg.Column([[sg.Text('Element 2')],[sg.Combo(elements[:nElements],default_value=elements[1],key='-el2-')]],vertical_alignment='t')]
            xLayout    = [sg.Column([[sg.Text('Start Element 2 Concentration')],[sg.Input(key='-xlo-',size=(inputSize,1))],
                          [sg.Text('Concentration unit')],[sg.Combo(['mole fraction'],default_value='mole fraction',key='-munit-')]],vertical_alignment='t'),
                          sg.Column([[sg.Text('End Element 2 Concentration')],[sg.Input(key='-xhi-',size=(inputSize,1))],
                          ],vertical_alignment='t'),
                          sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstep-',size=(8,1))]],vertical_alignment='t')]
            tempLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(inputSize,1))],
                          [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]],vertical_alignment='t'),
                          sg.Column([[sg.Text('End Temperature')],[sg.Input(key='-endtemperature-',size=(inputSize,1))],
                          ],vertical_alignment='t'),
                          sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstep-',size=(8,1))]],vertical_alignment='t')]
            presLayout = [sg.Column([[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(inputSize,1))],
                          [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]],vertical_alignment='t')
                          ]
            setupLayout = [elSelectLayout,xLayout,tempLayout,presLayout,[sg.Button('Run'), sg.Button('Refine', disabled = True), sg.Exit()]]
            setupWindow = sg.Window('Phase diagram setup', setupLayout, location = [400,0], finalize=True)
            while True:
                event, values = setupWindow.read(timeout=timeout)
                eventd, valuesd = dataWindow.read(timeout=timeout)
                if event == sg.WIN_CLOSED or event == 'Exit' or eventd == sg.WIN_CLOSED or eventd == 'Exit':
                    break
                elif event =='Run':
                    cancelRun = False
                    ntstep = values['-ntstep-']
                    nxstep = values['-nxstep-']
                    if (float(ntstep) * float(nxstep)) > 50000:
                        cancelRun = True
                        confirmLayout = [[sg.Text('The selected calculation is large and may take some time.')],[sg.Button('Continue'), sg.Button('Cancel')]]
                        confirmWindow = sg.Window('Large calculation confirmation', confirmLayout, location = [400,0], finalize=True)
                        while True:
                            event, values = confirmWindow.read(timeout=timeout)
                            if event == sg.WIN_CLOSED or event == 'Cancel':
                                break
                            elif event == 'Continue':
                                cancelRun = False
                                break
                        confirmWindow.close()
                    tlo= values['-temperature-']
                    thi = values['-endtemperature-']
                    pressure = values['-pressure-']
                    filename = 'inputs/pythonPhaseDiagramInput.ti'
                    tunit = values['-tunit-']
                    punit = values['-punit-']
                    xhi = values['-xhi-']
                    xlo = values['-xlo-']
                    el1 = values['-el1-']
                    el2 = values['-el2-']
                    if str(el1) == str(el2):
                        cancelRun = True
                        repeatLayout = [[sg.Text('The same element cannot be selected twice.')],[sg.Button('Cancel')]]
                        repeatWindow = sg.Window('Repeat element notification', repeatLayout, location = [400,0], finalize=True)
                        while True:
                            event, values = repeatWindow.read(timeout=timeout)
                            if event == sg.WIN_CLOSED or event == 'Cancel':
                                break
                            elif event == 'Continue':
                                cancelRun = False
                                break
                        repeatWindow.close()
                    munit = values['-munit-']
                    if cancelRun:
                        continue
                    writeInputFile(filename,xlo,xhi,tlo,thi,pressure,tunit,punit,munit,el1,el2,datafile)
                    ts = []
                    x1 = []
                    x2 = []
                    p1 = []
                    p2 = []
                    mint = 1e6
                    maxt = 0
                    mint, maxt = runCalc(ts, x1, x2, p1, p2, mint, maxt)
                    setupWindow.Element('Refine').Update(disabled = False)
                elif event =='Refine':
                    xRefLayout    = [sg.Column([[sg.Text('Start Concentration')],[sg.Input(key='-xlor-',size=(inputSize,1))]],vertical_alignment='t'),
                                  sg.Column([[sg.Text('End Concentration')],[sg.Input(key='-xhir-',size=(inputSize,1))]],vertical_alignment='t'),
                                  sg.Column([[sg.Text('# of steps')],[sg.Input(key='-nxstepr-',size=(8,1))]],vertical_alignment='t')]
                    tempRefLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperaturer-',size=(inputSize,1))]],vertical_alignment='t'),
                                  sg.Column([[sg.Text('End Temperature')],[sg.Input(key='-endtemperaturer-',size=(inputSize,1))]],vertical_alignment='t'),
                                  sg.Column([[sg.Text('# of steps',key='-tsteplabel-')],[sg.Input(key='-ntstepr-',size=(8,1))]],vertical_alignment='t')]
                    refineLayout = [xRefLayout,tempRefLayout,[sg.Button('Refine'), sg.Button('Cancel')]]
                    refineWindow = sg.Window('Phase diagram setup', refineLayout, location = [400,0], finalize=True)
                    while True:
                        event, values = refineWindow.read(timeout=timeout)
                        if event == sg.WIN_CLOSED or event == 'Cancel':
                            break
                        elif event =='Refine':
                            cancelRun = False
                            ntstep = values['-ntstepr-']
                            nxstep = values['-nxstepr-']
                            if (float(ntstep) * float(nxstep)) > 50000:
                                cancelRun = True
                                confirmLayout = [[sg.Text('The selected calculation is large and may take some time.')],[sg.Button('Continue'), sg.Button('Cancel')]]
                                confirmWindow = sg.Window('Large calculation confirmation', confirmLayout, location = [400,0], finalize=True)
                                while True:
                                    event, values = confirmWindow.read(timeout=timeout)
                                    if event == sg.WIN_CLOSED or event == 'Cancel':
                                        break
                                    elif event == 'Continue':
                                        cancelRun = False
                                        break
                                confirmWindow.close()
                            tlo = values['-temperaturer-']
                            thi = values['-endtemperaturer-']
                            xhi = values['-xhir-']
                            xlo = values['-xlor-']
                            if cancelRun:
                                continue
                            writeInputFile(filename,xlo,xhi,tlo,thi,pressure,tunit,punit,munit,el1,el2,datafile)
                            mint, maxt = runCalc(ts, x1, x2, p1, p2, mint, maxt)
            setupWindow.close()
        except:
            pass
