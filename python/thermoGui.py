import PySimpleGUI as sg
import subprocess
import math
import os
import sys

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

timeout = 50
inputSize = 20

windowList = []

class DataWindow:
    def __init__(self):
        windowList.append(self)
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
        self.folder = os.getcwd()+'/data'
        try:
            file_list = os.listdir(self.folder)
        except:
            file_list = []
        fnames = [
            f
            for f in file_list
            if os.path.isfile(os.path.join(self.folder, f))
            and f.lower().endswith((".dat", ".DAT"))
        ]
        fnames = sorted(fnames, key=str.lower)
        self.sgw = sg.Window('Thermochimica database selection', file_list_column, location = [0,0], finalize=True)
        self.sgw["-FILE LIST-"].update(fnames)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == "-FOLDER-":
            self.folder = values["-FOLDER-"]
            try:
                file_list = os.listdir(self.folder)
            except:
                file_list = []

            fnames = [
                f
                for f in file_list
                if os.path.isfile(os.path.join(self.folder, f))
                and f.lower().endswith((".dat", ".DAT"))
            ]
            fnames = sorted(fnames, key=str.lower)
            self.sgw["-FILE LIST-"].update(fnames)
        elif event == "-FILE LIST-":  # A file was chosen from the listbox
            try:
                datafile = os.path.join(
                    self.folder, values["-FILE LIST-"][0]
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
            except:
                return
            for el in elements:
                try:
                    index = atomic_number_map.index(el)+1 # get element indices in PT (i.e. # of protons)
                except ValueError:
                    if len(el) > 0:
                        if el[0] != 'e':
                            print(el+' not in list') # if the name is bogus (or e(phase)), discard
                    elements = list(filter(lambda a: a != el, elements))
            nElements = len(elements)
            if nElements == 0:
                return
            calcWindow = CalculationWindow(datafile,nElements,elements)
            self.children.append(calcWindow)

class CalculationWindow:
    def __init__(self, datafile, nElements, elements):
        windowList.append(self)
        self.datafile = datafile
        self.nElements = nElements
        self.elements = elements
        self.makeLayout()
        self.sgw = sg.Window(f'Thermochimica calculation: {os.path.basename(self.datafile)}', self.layout, location = [400,0], finalize=True)
        self.children = []
    def close(self):
        for child in self.children:
            child.close()
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()
        elif event == '-tdis-':
            self.sgw.Element('-endtemperature-').Update(visible = False)
            self.sgw.Element('-endtemperaturelabel-').Update(visible = False)
            self.sgw.Element('-ntstep-').Update(visible = False)
            self.sgw.Element('-tsteplabel-').Update(visible = False)
            self.sgw.Element('-pent-').Update(disabled = True)
            if values['-pent-']:
                self.sgw.Element('-pdis-').Update(value = True)
                self.sgw.Element('-endpressure-').Update(visible = False)
                self.sgw.Element('-endpressurelabel-').Update(visible = False)
        elif event == '-ten-':
            self.sgw.Element('-endtemperature-').Update(visible = True)
            self.sgw.Element('-endtemperaturelabel-').Update(visible = True)
            self.sgw.Element('-ntstep-').Update(visible = True)
            self.sgw.Element('-tsteplabel-').Update(visible = True)
            self.sgw.Element('-pent-').Update(disabled = False)
        elif event == '-pdis-':
            self.sgw.Element('-endpressure-').Update(visible = False)
            self.sgw.Element('-endpressurelabel-').Update(visible = False)
            self.sgw.Element('-pstep-').Update(visible = False)
            self.sgw.Element('-psteplabel-').Update(visible = False)
        elif event == '-pen-':
            self.sgw.Element('-endpressure-').Update(visible = True)
            self.sgw.Element('-endpressurelabel-').Update(visible = True)
            self.sgw.Element('-pstep-').Update(visible = True)
            self.sgw.Element('-psteplabel-').Update(visible = True)
        elif event == '-pent-':
            self.sgw.Element('-endpressure-').Update(visible = True)
            self.sgw.Element('-endpressurelabel-').Update(visible = True)
            self.sgw.Element('-pstep-').Update(visible = False)
            self.sgw.Element('-psteplabel-').Update(visible = False)
        elif event =='Run':
                temperature = values['-temperature-']
                pressure = values['-pressure-']
                filename = 'inputs/pythonInput.ti'
                masses = [0.0]*self.nElements
                for i in range(self.nElements):
                    if values['-'+self.elements[i]+'-'] != '':
                        masses[i] = values['-'+self.elements[i]+'-']
                tunit = values['-tunit-']
                punit = values['-punit-']
                munit = values['-munit-']
                if values['-ten-']:
                    tend = values['-endtemperature-']
                    ntstep = values['-ntstep-']
                else:
                    ntstep = 0
                if values['-pen-']:
                    pend = values['-endpressure-']
                    npstep = values['-pstep-']
                elif values['-pent-']:
                    pend = values['-endpressure-']
                    npstep = ntstep
                else:
                    npstep = 0
                with open(filename, 'w') as inputFile:
                    inputFile.write('! Python-generated input file for Thermochimica\n')
                    if float(ntstep) > 0:
                        tstep = (float(tend)-float(temperature))/float(ntstep)
                        inputFile.write('temperature          = ' + str(temperature) + ':' + str(tend) + ':' + str(tstep) + '\n')
                    else:
                        inputFile.write('temperature          = ' + str(temperature) + '\n')
                    if float(npstep) > 0:
                        pstep = (float(pend)-float(pressure))/float(npstep)
                        inputFile.write('pressure          = ' + str(pressure) + ':' + str(pend) + ':' + str(pstep) + '\n')
                    else:
                        inputFile.write('pressure          = ' + str(pressure) + '\n')
                    for i in range(self.nElements):
                        inputFile.write('mass(' + str(atomic_number_map.index(self.elements[i])+1) + ')           = ' + str(masses[i]) + '\n')
                    inputFile.write('temperature unit         = ' + tunit + '\n')
                    inputFile.write('pressure unit          = ' + punit + '\n')
                    if values['-pent-']:
                        inputFile.write('step together     = .TRUE.\n')
                    inputFile.write('mass unit         = ' + munit + '\n')
                    inputFile.write('data file         = ' + self.datafile + '\n')
                    inputFile.write('print mode        = 2\n')
                    inputFile.write('debug mode        = .FALSE.\n')
                    if values['-json-']:
                        inputFile.write('write json     = .TRUE.\n')
                thermoOut = subprocess.check_output(['./bin/ThermochimicaInputScriptMode',filename]).decode("utf-8")
                nLines = thermoOut.count('\n')
                if (nLines < 5000):
                    resultOutput = [[sg.Column([[sg.Multiline(thermoOut, size = (65, nLines))]], size = (400, 800), scrollable = True, vertical_scroll_only = True)]]
                else:
                    resultOutput = [[sg.Text('Output is too large to display')]]
                resultWindow = ResultWindow(resultOutput)
                self.children.append(resultWindow)
    def makeLayout(self):
        tempLayout = [sg.Column([[sg.Text('Temperature')],[sg.Input(key='-temperature-',size=(inputSize,1))],
                      [sg.Text('Temperature unit')],[sg.Combo(['K', 'C', 'F'],default_value='K',key='-tunit-')]],vertical_alignment='t'),
                      sg.Column([[sg.Text('End Temperature',key='-endtemperaturelabel-',visible=False)],
                      [sg.Input(key='-endtemperature-',size=(inputSize,1),visible=False)],
                      [sg.Text('Temperature range:',visible=False)],
                      [sg.Radio('Disabled', 'trange', default=True, enable_events=True, key='-tdis-')],
                      [sg.Radio('Enabled', 'trange', default=False, enable_events=True, key='-ten-')]],vertical_alignment='t'),
                      sg.Column([[sg.Text('# of steps',key='-tsteplabel-',visible=False)],[sg.Input(key='-ntstep-',size=(8,1),visible=False)]],vertical_alignment='t')]
        presLayout = [sg.Column([[sg.Text('Pressure')],[sg.Input(key='-pressure-',size=(inputSize,1))],
                      [sg.Text('Pressure unit')],[sg.Combo(['atm', 'Pa', 'bar'],default_value='atm',key='-punit-')]],vertical_alignment='t'),
                      sg.Column([[sg.Text('End Pressure',key='-endpressurelabel-',visible=False)],[sg.Input(key='-endpressure-',size=(inputSize,1),visible=False)],
                      [sg.Text('Pressure range:')],
                      [sg.Radio('Disabled', 'prange', default=True, enable_events=True, key='-pdis-')],
                      [sg.Radio('Enabled', 'prange', default=False, enable_events=True, key='-pen-')],
                      [sg.Radio('Enabled, step\nwith temperature', 'prange', default=False, enable_events=True, key='-pent-')]],vertical_alignment='t'),
                      sg.Column([[sg.Text('# of steps',key='-psteplabel-',visible=False)],[sg.Input(key='-pstep-',size=(8,1),visible=False)]],vertical_alignment='t')
                      ]
        elem1Layout = [[sg.Text('Composition 1')]]
        elem2Layout = [[sg.Text('Composition 2')]]
        for el in self.elements:
            elem1Layout.append([sg.Text(el)])
            elem1Layout.append([sg.Input(key=f'-{el}1-',size=(inputSize,1))])
        for el in self.elements:
            elem2Layout.append([sg.Text(el)])
            elem2Layout.append([sg.Input(key=f'-{el}2-',size=(inputSize,1))])
        if (self.nElements < 8):
            self.layout = [tempLayout,
                          presLayout,
                          [sg.Column(elem1Layout),sg.Column(elem2Layout,key='composition2',visible=False)],
                          [sg.Text('Mass unit')],
                          [sg.Combo(['moles', 'kg', 'atoms', 'g'],default_value='moles',key='-munit-')],
                          [sg.Checkbox('Save JSON',key='-json-')],
                          [sg.Button('Run'), sg.Exit()]]
        else:
            self.layout = [tempLayout,
                          presLayout,
                          [sg.Column(elem1Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True),
                           sg.Column(elem2Layout,vertical_alignment='t', scrollable = True, vertical_scroll_only = True, expand_y = True,key='composition2',visible=False)],
                          [sg.Text('Mass unit')],
                          [sg.Combo(['moles', 'kg', 'atoms', 'g'],default_value='moles',key='-munit-')],
                          [sg.Checkbox('Save JSON',key='-json-')],
                          [sg.Button('Run'), sg.Exit()]]

class ResultWindow:
    def __init__(self, layout):
        windowList.append(self)
        self.sgw = sg.Window('Thermochimica output',layout, location = [825,0], finalize=True)
    def close(self):
        self.sgw.close()
        if self in windowList:
            windowList.remove(self)
    def read(self):
        event, values = self.sgw.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            self.close()


if not(os.path.isfile('bin/ThermochimicaInputScriptMode')):
    errorLayout = [[sg.Text('No Thermochimica executable available.')],
                   [sg.Text('Either Thermochimica has not been built (run make),')],
                   [sg.Text('or this script was not executed from Thermochimica root directory.')],
                   [sg.Button('Exit')]]
    errorWindow = sg.Window('Thermochimica Error Message', errorLayout, location = [0,0], finalize=True)
    while True:
        event, values = errorWindow.read(timeout=timeout)
        if event == sg.WIN_CLOSED or event == 'Exit':
            break
    errorWindow.close()
    sys.exit()

dataWindow = DataWindow()
while len(windowList) > 0:
    for window in windowList:
        window.read()
