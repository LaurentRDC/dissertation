Note: this script was modified from https://gitlab.com/VADemon/onedeath for use in Github Actions.

## OneDeath - a crappy OneDrive shared folder downloader

OneDeath will download entire public OneDrive/Sharepoint folders (shares) as a single batch, using wget. It will check SHA1-hashes after downloading remote to ensure correctness. (doesn't work for business drives)

I wrote this to overcome onedrive's bulk downloading """capabilities""" It just uses a bunch of hastily written Lua code to make wget do the right things with the API.

If you have any problems, open an issue and ask.

### Usage:

`lua main.lua '<SHARE URL>'` - use single or double quotes around the URL!

**Accepted URLs:**

Only use the shortened URLs! Don't copy the URL from an opened browser tab: they don't grant access tokens!

* OneDrive: `https://1drv.ms/f/s!<longText>`

* Sharepoint: `https://<site>.sharepoint.com/:f:/g/personal/<drive>/<longText>` - not sure about :f: and personal

Optional argument only accepted as argument `#2`: `lua main.lua '<SHARE URL> [-p|-w]`

**Additional options:**

* Only print wget links, no DL: `-w`
    * extracting should *currently* work with `<command> | grep -C 1 'EXTRACTED DOWNLOAD URL' | grep ^wget`

* [Debug] only convert and print the short `1drv.ms` URL to an API URL: `-p`

**Want to stop the script and abort download?**

*Ctrl+C* won't work because it's only killing the foreground task :) 

1. Use *Ctrl+Z* to halt the script, it will display a *[number]* on the left (or run `jobs` to see all)
2. Run `kill %<number>` to terminate it (e.g. `kill %1`)

#### Example output:

```
(6/8) File:  15.0 GB Somefile.zip
Downloading to: <SharerName>_<ShareFolder>/<File>.zip
Local file [ 14.5 GB] is smaller than remote [ 15.0 GB] (96.6%), may be able to resume download
<File>.zip  96%[++++++++++++++++++++++++++++++++++++++++++    ]  14.49G  2.78MB/s    eta 3m 11s
```

### Requirements:

* Lua
* wget
* grep
* coreutils: `sha1sum`, `touch`, `stat`
   * [Line 3:](https://gitlab.com/VADemon/onedeath/blob/master/main.lua#L3) Change `checkSha1` to false if you wish to disable hash checks

### Installation:

#### Linux/Debian

1. `sudo apt-get install -y wget lua`

#### Windows

* Get [Cygwin](https://www.cygwin.com/) and install `wget` & `lua`, `coreutils` packages there

* Then to launch the downloader, open Cygwin Terminal, navigate to the folder with `main.lua`

