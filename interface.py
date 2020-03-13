
from tkinter import filedialog
from tkinter import messagebox
from tkinter import *
import os
import re
import helpText
import webbrowser



def initGui():
    cwd = os.getcwd()

    #initialize main window
    root = Tk() 
    root.title("ComplexBuilder") 
    img = Image("photo", file='CB_icon.png')
    root.call('wm','iconphoto',root._w, img)
    #root.geometry("400x400")
    #root.option_add("*Button.Background", "grey")
    #root.option_add("*Button.Foreground", "#4f91b4")
    root.configure(background = "#2da6e6")
    
    
    def simpleName(inputPath):
        """Return the last element of a path to display just the input name"""
        justInput=re.search('([^/]+)/?$', inputPath)   
        return justInput.group(0)  
    
    def fastaSelection():
        """Open window with *.fa files in the current directory to select the FASTA file"""
        filetext="Select"
        root.filename =  filedialog.askopenfilename(initialdir = cwd,title = "Select the FASTA file",
                                                    filetypes = (("FASTA files","*.fa"),
                                                                ("all files","*.*")))
        print (root.filename)
        bF["text"]= simpleName(root.filename) if root.filename else filetext
        return root.filename

    def pdbSelection():
        """Open window with directories present in the current folder to select the the one 
        with the pair interactions in PDB"""
        filetext="Select"
        root.directory = filedialog.askdirectory(title = "Select the directory with PDB files")
        print (root.directory)
        bP["text"]= simpleName(root.directory) if root.directory else filetext
        return root.directory 

    def oSelection():
        """Store the name introduced for the output file """
        filetext="Confirm"
        print (outputName.get())
        bO["text"]= simpleName(outputName.get()) if outputName.get() else filetext
        return outputName.get()

    def stSelection():
        """Open window with files in the current directory to select the the stoichiometry file"""
        filetext="Select"
        root.filenameSt =  filedialog.askopenfilename(initialdir = cwd,title = "Select the stoichiometry file",
                                                    filetypes = (("Text files","*.txt"),
                                                                ("all files","*.*")))
        print (root.filenameSt)
        frameSt.pack(pady="15", padx="50",anchor=W, fill=X, expand=YES)
        bSt["text"]= simpleName(root.filenameSt) if root.filenameSt else filetext   
        labelO["text"]="Output directory name:                       "
        labelSt.lift()
        labelSt.pack(side=LEFT )
        bSt.lift()
        bSt.pack()
        return root.filenameSt

    def checkVerbose():
        """return TRUE if the log checkbox if checked"""
        return verb.get()

    def showHelp():
        """Display help of the ComplexBuilder in a new window"""
        t = Toplevel(root)
        t.wm_title("Help")
        l = Label(t, text=helpText.helpMessage)
        l.pack(side="top", fill="both", expand=True, padx=100, pady=100)

    def openTutorial():
        """Display README file from gitHub"""
        webbrowser.open_new("https://github.com/Argonvi/SBI-PYT_Project/blob/master/README.md")

    def runCB():
        """Check that all the mandatory inputs are defined and run the program"""
        if outputName.get()!='':
            try:
                print(root.filename)
                print(root.directory )
                root.destroy()
            except:
                messagebox.showerror("Ups!", "You should select the input FASTA file and the directory with the PDB interaction files.")
                
        else:
            messagebox.showerror("Ups!", "You should specify and confirm the name of the directory where the output files will be stored.")
             


    menubar = Menu(root, tearoff=0)

    filemenu = Menu(menubar)
    filemenu.add_command(label="Add stoichiometry", command=stSelection)
    #filemenu.add_command(label="Energy analysis", command=menu_action)
    #filemenu.add_command(label="Open lof file", command=menu_action)
    filemenu.add_separator()
    filemenu.add_command(label="Quit", command=root.quit)

    helpmenu = Menu(menubar)
    helpmenu.add_command(label="Help", command=showHelp)
    helpmenu.add_separator()
    helpmenu.add_command(label="About", command=openTutorial)

    menubar.add_cascade(label="Options", menu=filemenu)
    menubar.add_cascade(label="Help", menu=helpmenu)


    root.config(menu=menubar)



    #Console frame
    frame = Frame(root)
    frame.pack(pady="10", padx="10")
    
    wellcome = Label(frame, text="Introduce the FASTA file and the directory with the interactions in PDB:")
    wellcome.pack(side=TOP,pady="15", padx="45")
    

    # FASTA widget
    frameFA = Frame(frame)
    frameFA.pack(pady="15", padx="50",side=TOP,anchor=W, fill=X, expand=YES)
    labelF = Label(frameFA, text="Sequence file:                                    ")
    labelF.pack( side=LEFT )
    #Button FASTA file
    bF = Button(frameFA, text="Select", command=fastaSelection)
    bF.pack( )


    #PDB widget
    framePDB = Frame(frame)
    framePDB.pack(pady="15", padx="50",anchor=W, fill=X, expand=YES)
    labelF = Label(framePDB, text="Interactions directory:                       ")
    labelF.pack( side=LEFT )
    #Button PDB dir
    bP = Button(framePDB, text="Select", command=pdbSelection)
    bP.pack(  )


    #Output widget
    frameO = Frame(frame)
    frameO.pack(pady="15", padx="50",anchor=W, fill=X, expand=YES)
    labelO = Label(frameO, text="Output directory name:    ")
    labelO.pack( side=LEFT )
    #textBox stoichiometry
    outputName = StringVar()
    nameEntered = Entry(frameO, width = 15, textvariable = outputName)
    nameEntered.pack(side=LEFT, padx="5")
    bO = Button(frameO, text = "Confirm", command = oSelection)
    bO.pack(side=LEFT)        


    #Verbose widget
    frameV = Frame(frame)
    frameV.pack(pady="15", padx="50",side=LEFT,anchor=W, fill=X, expand=YES)
    labelF = Label(frameV, text="Do you want to create a log file?")
    labelF.pack( side=LEFT )
    #Checkbox verbose
    verb = BooleanVar()
    cV = Checkbutton(frameV, text="I do!", variable=verb, command=checkVerbose)
    cV.pack()

    #stoichiometry widget
    frameSt = Frame(frame)
    #frameSt.pack(pady="15", padx="50",anchor=W, fill=X, expand=YES)
    frameSt.pack_forget()
    labelSt = Label(frameSt, text="Define a stoichiometry:  ")
    #labelSt.pack( side=LEFT )
    labelSt.lower()
    bSt = Button(frameSt, text = "Select", command = stSelection)
    #bSt.pack()
    bSt.lower()
    


    #Submit widget
    bEnter = Button(root, text="Build complex", command=runCB, bd="2", 
                     relief="solid",padx="10", pady="10") #fg="red", bg="blue",     
    bEnter.pack( side=BOTTOM ,pady="10", padx="1") 
    root.mainloop()




    #Parse inputs 
    inputs=[]
    inputs.append( simpleName(root.filename)) # FASTA
    inputs.append( simpleName(root.directory)) # PDB dir
    inputs.append(checkVerbose()) # verbose T/F
    try: # stoichiometry None/file
        inputs.append( simpleName(root.filenameSt))
    except:
        inputs.append(None)
    inputs.append(outputName.get())
    print(inputs)
    return inputs

    

