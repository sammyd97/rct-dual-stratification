# How to Push to GitHub

Your repository is ready! Follow these steps to push to GitHub:

## Step 1: Create a new repository on GitHub

1. Go to https://github.com/new
2. Repository name: `rct-dual-stratification`
3. Description: `RCT facility randomization tool with dual stratification (LGA + Volume)`
4. Choose Public or Private
5. **DO NOT** initialize with README, .gitignore, or license (we already have these)
6. Click "Create repository"

## Step 2: Push your code

After creating the repository, run these commands in Terminal:

```bash
cd "/Users/skhsa/Desktop/RCT randomisation"

# Add your GitHub repository as remote
git remote add origin https://github.com/sammyd97/rct-dual-stratification.git

# Push to GitHub
git branch -M main
git push -u origin main
```

## Alternative: If you have SSH set up

```bash
git remote add origin git@github.com:sammyd97/rct-dual-stratification.git
git branch -M main
git push -u origin main
```

## What's included in this repository:

✅ `randomize_facilities.py` - Main randomization script (supports dual stratification)  
✅ `README.md` - Complete documentation  
✅ `Benue_randomised.csv` - Example dual-stratified results (LGA + Volume)  
✅ `Benue_balance_report.csv` - Example balance report  
✅ `.gitignore` - Git ignore rules  

Ready to use with your own data!

---

Need help? Just create the repository on GitHub and follow the commands above!

