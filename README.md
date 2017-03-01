# CX-4220-PA-1
Programming Assignment 1

Setting up git:
Download git via command line tools

$ git config --global user.name "YOUR NAME"
$ git config --global user.email "YOUR EMAIL USED FOR GITHUB"

Navigate to where you want to clone the repository (where you want the folder containing program files to go).
This part might be weird

$ git clone https://github.gatech.edu/ttippens6/CX-4220-PA-1.git

You now have a local copy you can work on.

Using git:
Before making changes, do

$ git pull
$ git checkout master

Then do your changes.
Once your changes have been added, you 'commit them'. Commit once, after your done working.
Do:

$ git commit -a -m "Brief description of the changes you made"

After the commit, you'll need to push your changes to the online repo:

$ git push

and you're done.
