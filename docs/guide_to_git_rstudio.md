# Git guide

Here's how to get started with git integration with RStudio and basics of Rstudio.

### Installing Git and interfacing it with RStudio

Here are some instructions provided by the folks at RStudio.

https://support.rstudio.com/hc/en-us/articles/200532077-Version-Control-with-Git-and-SVN

### Accessing the git repo in RStudio

In RStudio, on the top menu, select "File > New Project". You should get a box with three options - "New Directory", "Existing Directory", or "Version Control". Select "Version Control". From there, select "Git". You should get a dialog box with three field - "Repository URL", "Project directory name", and "Create as a subdirectory of".

Now, navigate to the repo page on Github. On the right-hand side of the page, there's a button for "Clone or download". *Cloning* a repository puts the repository (i.e., all of the files and folders) on your computer in its current state. If you click on this button, it will give you a URL. Copy this URL.

Now, in the dialog box in RStudio, copy and paste this URL into the "Repository URL" field. This URL should end in ".git". The "project directory name" field should autopopulate with the name of the repository. The "project directory name" is just the name of the folder in which your repository will be stored. You can rename the repository whatever you'd like, but RStudio will suggest the name from the git file. The final field is the directory on your computer where you would like the repository to be downloaded and stored.

Once you click the "create project" button, it will clone the repository to the directory you specified. You'll note that there will now be an icon for this project in the upper-right corner of your RStudio window with the name of the repository you specified. You'll also have a new panel window in RStudio for "Git".

### Using features of Git

The idea behind git is it provides a useful and organized way for multiple people to simultaneously work on a project. The most recent version of the repository is always stored on GitHub. The version of the repository stored on your computer is a **local** version of the repository. You can sync up your local repository to match the one online by *pulling* the online repository into your repo. Likewise, you can add features or changes in your repository to the online version of the repo (so that other people can use those features) by **pushing** your version of the repository into the online one.

##### Step 1: Pulling in RStudio

If there are changes to the current repository you want to add to your repository (e.g., data files added or scripts updated), then you can add them to your local repository by pulling.

Do this in RStudio by navigating to the "Git" panel in RStudio. At the top of the panel there are a series of buttons, (e.g., "Diff", "Commit"). *Pull by clicking on the downward-pointing arrow.* This is the only step required for pulling.

Pulling re-downloads new or updated files from the entire online repository to your system. If you have locally stored files which are not on the repo, *they will not be overwritten*. So, you can have work in progress and have it be unaltered before each pull. However, if you have work in a file *that already exists in the repository* (or a file with a name of a file already present in the repo), then *pulling will wipe out your changes if they aren't committed*. (Details on what a commit is below).

One final note is that it's a good idea to **pull often**. It's good practice to pull every time you're going to start working on something in the repository. It's also a good idea to **pull before every time you push**.

##### Step 2: Committing in RStudio

The first step to pushing things onto the repository is storing changes as *commits*. Git stores individual changes as commits, stores these commits in what's called a "staging area". Then you push your commits all at once to the online repository.

Once you have a change made to a file, the file should appear in the "Git" panel of RStudio. A brand new file should have two yellow squares next to it, while a pre-existing (in the repo) file will have a blue square next to it. To commit your changes, click on the white button next to the name of the file, then click on the "Commit" button. This will open a box where all of the things you have added to the file will be in green, while all of the things removed from the file are in red. 

With every commit, you should add a *commit message*. Every commit is stored in the repository (this allows you to go back to previous versions of the repository online) - commit messages are short messages explaining what is happening in each commit. This is very useful if you need to go back to restore a previous version of the repository. An example of a commit message is "adding code for cleaning data", or "initializing script to draw figures".

Once you have written a good commit message and the file is "staged" (i.e., the white "stage button" for that file is filled in), press the "Commit" button.

You can change multiple files with one commit by clicking the "staged" button for several files. You can also put together multiple commits at a time before pushing.

Every GitHub repo puts the list of prior commits online. Here is the list of commits for this repo: https://github.com/EBIO6100Spring2020/Data-sandbox/commits/master. You can find this list by going to the main page for the repository and hitting the "Commits" button below the name of the repository.

##### Step 3: Pushing in RStudio

Once you have at least one commit ready to put into the online repo, all you need to do to commit is hit the "Push" button in RStudio. This is the green upward-facing arrow. **Remember to pull before you push.** This will now add your committed changes to the online (master) version of the repo. Note that changes you made that you did not stage and commit will not be reflected in the repo.


