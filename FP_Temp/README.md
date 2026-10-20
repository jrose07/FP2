###### How to Use Template ######

1. Clone Repository and rename it to current Experiment
2. Use "git clone <https of new Repo VXX>" to save new Repo into local files
3. Use "make NAME='Experiment Name' " to rename inside the files
4.  Because of Some weird thing you then have to type 'rm ./.git/index' in order to commit the rename thing.
5. Write the Python Script inside of "plot.py" 
6. Write the Sections inside "vXXX/content"

###### Use make in different Modes #####

-"make MODE=write" lets you write inside the latex files while the pdf is simultaneously updating.
-"make plots" will only run the python file
-"make" will run everything 


Hope that was everything ;)
