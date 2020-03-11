
from tkinter import filedialog
from tkinter import *
import os
import re



def initGui():
    cwd = os.getcwd()
    root = Tk() #initialize window
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

    def stSelection():
        """Store the stoichiometry entered by the user"""
        labelSt.configure(text= 'Stoichiometry: ' + st.get())

    def checkVerbose():
        """return TRUE if the log checkbox if checked"""
        return verb.get()






    #Console frame
    frame = Frame(root)
    frame.pack(pady="10", padx="10")
    
    wellcome = Label(frame, text="Introduce the FASTA file and the directory with the interactions in PDB:")
    wellcome.pack(side=TOP,pady="15", padx="45")
    

    # FASTA widget
    frameFA = Frame(frame)
    frameFA.pack(pady="15", padx="50",side=TOP,anchor=W, fill=X, expand=YES)
    labelF = Label(frameFA, text="Sequence file:            ")
    labelF.pack( side=LEFT )
    #Button FASTA file
    bF = Button(frameFA, text="Select", command=fastaSelection)
    bF.pack( )


    #PDB widget
    framePDB = Frame(frame)
    framePDB.pack(pady="15", padx="50",anchor=W, fill=X, expand=YES)
    labelF = Label(framePDB, text="Interactions directory:")
    labelF.pack( side=LEFT )
    #Button PDB dir
    bP = Button(framePDB, text="Select", command=pdbSelection)
    bP.pack(  )


    #stoichiometry widget
    frameSt = Frame(frame)
    frameSt.pack(pady="15", padx="50",anchor=W, fill=X, expand=YES)
    labelSt = Label(frameSt, text="Define a stoichiometry (optional):")
    labelSt.pack( side=LEFT )
    #textBox stoichiometry
    st = StringVar()
    button = Button(frameSt, text = "Apply", command = stSelection)
    button.pack(pady="15", padx="17",side=RIGHT)
    stEntered = Entry(frameSt, width = 10, textvariable = st)
    stEntered.pack(pady="15", padx="5",side=RIGHT)
 
    


    #Verbose widget
    frameV = Frame(frame)
    frameV.pack(pady="15", padx="50",side=LEFT,anchor=W, fill=X, expand=YES)
    labelF = Label(frameV, text="Do you want to create a log file?")
    labelF.pack( side=LEFT )
    #Checkbox verbose
    verb = BooleanVar()
    cV = Checkbutton(frameV, text="I do!", variable=verb, command=checkVerbose)
    cV.pack(  )
    


    #Submit widget
    bEnter = Button(root, text="Build complex", command=root.destroy, bd="2", 
                     relief="solid",padx="10", pady="10") #fg="red", bg="blue",     
    bEnter.pack( side=BOTTOM ,pady="10", padx="1") 
    root.mainloop()


    #Parse inputs 
    inputs=[]
    inputs.append( simpleName(root.filename)) 
    inputs.append( simpleName(root.directory))
    inputs.append(checkVerbose())
    print(inputs)
    return inputs

    

