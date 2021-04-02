#!/usr/bin/env lua
---
checkSha1 = true
showDownloadProgress = false -- enables wget progress bar (only for files)
MAX_PATH = 130 -- shorten some folders/names if they're over this limit... not everyone has a good FS with > 255 char limit
---

-- Sharepoint access tokens only valid for 1h: 
-- "Sorry, something went wrong / The access token has expired. It's valid from '1/21/2020 12:39:02 AM' and to '1/21/2020 1:39:02 AM'."

require("lib.base64")
JSON = require("lib.JSON")

userAgent = "Mozilla/5.0 (X11; Linux x86_64; rv:71.0; Microsoft, I hate you111) Gecko/20100101 Firefox/71.0"
cookiesFile = "cookies.txt"

cookiesOptions = "--keep-session-cookies --save-cookies ".. cookiesFile .." --load-cookies ".. cookiesFile
wgetOptions = cookiesOptions
			.. " --continue " -- dont restart download, try to continue. Absolutely required for Sharepoint large files
			.. " --retry-connrefused --waitretry=1 --read-timeout=20 --timeout=15 -t 10 "
			.. " -U '".. userAgent .."' "
			--.. " --header='Accept: application/json;odata=verbose' --header='Content-Type: application/json;odata=verbose'"

sizeSuffix = {"B", "KB", "MB", "GB", "TB", "PB", "EB", "ZB", "YB"}
function formatSize(s, lpadding)
	lpadding = lpadding or 1
	suffixCnt = 1
	while s > 999 and suffixCnt < 10 do
		s = s / 1024
		suffixCnt = suffixCnt + 1
	end
	return string.format("% ".. lpadding ..".1f %s", s, sizeSuffix[suffixCnt])
end

do
	local magicChars = [=[.-()?%*[]+$^]=]
	local magicPattern = "[".. magicChars:gsub(".", "%%%1") .."]"

	string.regexEscape = function (str)
		return (str:gsub(magicPattern, "%%%1")) -- only first argument returned
	end
end

-- all escapes in bash paths (includes space)
--!"$&'( )*,:;<=>?@[\]^`{|}
function string.shellEscapeLiteral(str) -- escapes string that will be put inside single quotes like 'I kissed john's hamster'
	return (str:gsub("'", "'\\''"))
end

function string.shellMakeLiteral(str) -- escapes string and surrounds in ' single-quotes - ready to be put into shell calls
	return "'".. str:shellEscapeLiteral() .."'"
end

function sanitisePathLength(str, fromRight)
	local newStr = str
	
	for name in str:gmatch("[^./\\]+") do -- dont break paths
		if #name > MAX_PATH then
			if fromRight then
				--print(name:sub(1, MAX_PATH):gsub("%%","%%%%"))
				newStr = newStr:gsub(name:regexEscape(), 
					name:sub(1, MAX_PATH):gsub("%%","%%%%") .. "CUT", 1)
			else
				newStr = newStr:gsub(name:regexEscape(), "CUT" .. name:sub(-MAX_PATH):gsub("%","%%"), 1)  -- return end of string
			end
		end
	end
	
	return newStr
end

-- Lua5.2+ os.execute-style command
if _VERSION == "Lua 5.1" or _VERSION == "Lua 5.0" then
	local execute = os.execute
	function os.execute(...)
		local status = execute(...)
		if status == 0 then
			return true, "exit", 0
		else
			-- the values were shifted by <<8 bits prior to 5.2 release...
			return nil, "exit", status and status/256 or status
		end
	end
end

function mkdir(path)
	os.execute("mkdir -p ".. path:shellMakeLiteral())
end

function sleepRetry(retries)
	local sleepTime = math.ceil(retries^1.5 * 5 + 5)
	print("Retrying in ".. sleepTime .."...")
	os.execute("sleep ".. sleepTime)
end

function wget(url)
	url = url:gsub("'", "%%27")
	return wgetEx(url, "-O -")
end

function wgetEx(url, extraOptions)
	local p = io.popen("wget ".. wgetOptions .." -nv ".. url:shellMakeLiteral() .." ".. extraOptions, "r")
	local s = p:read("*a")
	p:close()
	return s
end

function getRedirect(url)
	return wgetEx(url, "--max-redirect=0 -S -O /dev/null 2>&1"):match("Location: ([%S]+)")
end

-- returns diff wall time needed for dl
skipDl = false
function wgetFile(url,toPath, extraOptions, retryLimit)
	extraOptions = extraOptions or ""
	if showDownloadProgress then
		extraOptions = extraOptions .. " --progress=bar:force:noscroll --show-progress "
	end
	
	local time, exitcode = wgetFileEx(url, toPath, "-q -nv ".. extraOptions)
	
	local retries, retryMax = 0, math.max(1, retryLimit or 10)
	while retries < retryMax do
		retries = retries + 1
		
		if exitcode ~= 0 then
			print("Non zero exit code (".. exitcode ..") returned:\nFor Url: ".. url .." and targetPath: ".. toPath)
			
			if exitcode == 4 then -- server-closed connection, maybe finished
				print("Server closed connection, not initiating inner retry (probably finished).")
				return false, time
			elseif exitcode == 6 then -- 401 Unauthorized
				print(retryMax.."retryLimit set as:", retryLimit)
				return false, time
			end
		elseif exitcode == 0 then
			break
		end
		
		sleepRetry(retries)
		time, exitcode = wgetFileEx(url, toPath, "-v")
	end
	if retries == retryMax then
		--error("Max retry count reached! Stopping")
		return false, time
	end
	
	return true, time
end

function wgetFileEx(url,toPath, extraOptions)
	local start = os.time()
	url = url:gsub("'", "%%27")
	
	-- catchsegv ./wget
	local _,_,exitcode = os.execute("wget ".. extraOptions .." ".. wgetOptions .." -O "..toPath:shellMakeLiteral().." ".. url:shellMakeLiteral())
	local diff = os.time()-start
	return diff~=0 and diff or 1, exitcode
end

function genSha1(path)
	local p
	repeat
		p = io.popen("sha1sum --binary ".. path:shellMakeLiteral(), "r")
		if not p then print("sha1 pipe doesnt exist!") end
	until p

	local s = p:read("*a")
	p:close()
	return s:sub(1,40)
end

-- extracts cookies from Netscape-type cookie file
function extractCookies(path)
	local file = io.open(path, "rb")
	local cookies = {
		-- [name] = "cookieval"
	}
	for line in file:lines() do
		if #line > 1 and not line:find("^[ \t\v]*#") then
			--<site>.sharepoint.com	FALSE	/	TRUE	0	FedAuth	77u/PD94bWwgd
			site, _, _, _, _, name, value = line:match("^([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)\t([^\t]+)$")
			if not name or not value then
				print("Couldnt extract cookie from file!")
				print(line)
			else
				cookies[name] = value
			end
		end
	end
	file:close()
	return cookies
end

-- returns wget-compatible args for cookies headers from Netscape-type cookie file
function genCookieHeaders(path)
	local cookies = extractCookies(path)
	local headerArgs = {}
	for name, val in pairs(cookies) do
		table.insert(headerArgs, string.format(" --header 'Cookie: %s=%s'", name, val))
	end
	return table.concat(headerArgs, " ")
end

function getFileSize(path)
	local p = io.popen("stat -c '%s' ".. path:shellMakeLiteral(), "r")
	local text = p:read("*l")
	p:close()
	
	local size = tonumber(text) or 0
	return size
end

function getFileMtime(path)
	local p = io.popen("stat -c '%Y' ".. path:shellMakeLiteral(), "r")
	local text = p:read("*l")
	p:close()
	
	local mtime = tonumber(text)
	return mtime
end

function setFileTime(path, dateStr)
	os.execute("touch --time=modify -d ".. dateStr:shellMakeLiteral() .." ".. path:shellMakeLiteral())
end

-- returns epoch time
function dateTzToTime(str)
	-- "2111-45-12T12:00:03.123Z"
	local p = io.popen("date --date ".. str:shellMakeLiteral() .." +%s", "r")
	local text = p:read("*l")
	p:close()
	return tonumber(text)
end

function getFolderName(apiUrl)
	local rootJson = wget(apiUrl:gsub("(/root/).+$", "%1"))	
	local tbl = JSON:decode(rootJson)

	if tbl.name and tbl.createdBy and tbl.createdBy.user and tbl.createdBy.user.displayName then
		return tbl.createdBy.user.displayName:gsub(" ", "") .. "_" .. tbl.name
	else
		return apiUrl:gsub(".+!",""):gsub("[^A-Za-z0-9_%-]", "_"):sub(-32)
	end
end

function fileExists(path)
	local f = io.open(path, "r")
	if f then
		f:close()
		return true
	else
		return false
	end
end

function printWgetCommand(fileName, path, dlpath, url)
	local cookieArg = genCookieHeaders(cookiesFile)
	print("File: ".. fileName)
	print("EXTRACTED DOWNLOAD URL:")
	print(("wget %s %s -O %s %s"):format(wgetOptions, cookieArg, 
			dlpath:sub(1+select("2",dlpath:find(path,1,true))):shellMakeLiteral(),
			url:shellMakeLiteral()
		)
	)
end

-- returns true = success
-- false = failed
function apiDownloadFile(item, path, dlpath, wgetRetries)
	local dltime, ok = -1
	local downloadUrl = sharepointDownload and (item["@odata.id"].."/content") or item["@content.downloadUrl"]
	if extractLinksOnly then
		printWgetCommand(item.name, path, dlpath, downloadUrl)
		return true
	end
	
	local downloadNeeded = true
	local hashingNeeded = true
	print("Downloading to: ".. dlpath)
	
	local localSize
	if fileExists(dlpath) then
		localSize = getFileSize(dlpath)
		if localSize < item.size then
			print(string.format("Local file [%s] is smaller than remote [%s] (%04.1f%%), may be able to resume download",
			formatSize(localSize,5), formatSize(item.size,5), localSize/item.size*100))
			
		elseif localSize == item.size then
			print("Local file exists and size is correct!")
			
			local localMtime = getFileMtime(dlpath)
			local remoteMtime = dateTzToTime(item.lastModifiedDateTime)
			local localDate = os.date("!%Y-%m-%dT%H:%M:%SZ", localMtime)
			local remoteDate = os.date("!%Y-%m-%dT%H:%M:%SZ", remoteMtime)
			
			if localMtime >= remoteMtime then
				downloadNeeded = false
				
				if localMtime > remoteMtime then
					print("Local date is newer than remote date! This shouldn't happen, check your files:")
					print(string.format("(UTC) Local: %d %s, Remote: %d %s", localMtime, localDate, remoteMtime, remoteDate))
					--setFileTime(dlpath, item.lastModifiedDateTime) -- force-update mismatching dates
				end
			else
				print("Remote file was updated! (UTC) Local date: ".. localDate, "Remote date: ".. remoteDate)
			end
		end
	end
	local oldSize = localSize or 0
	
	if downloadNeeded then
		ok, dltime = wgetFile(downloadUrl, dlpath, wgetRetries)
		localSize = getFileSize(dlpath)
		if not ok and localSize ~= item.size then
			print("Not yet finished, took ~".. dltime .."s, ".. formatSize((localSize-oldSize)/dltime,5) .."/s")
			return false
		end
	end
	
	local localSize = localSize or getFileSize(dlpath)
	if localSize ~= item.size then
		print("Size mismatch! Expected ".. item.size ..", but local is ".. localSize)
		if localSize > item.size then os.exit(1) end
		return false
	end
	
	if dltime > -1 then
		setFileTime(dlpath, item.lastModifiedDateTime)
		io.stdout:write("Done in ~".. dltime .."s, ".. formatSize((localSize-oldSize)/dltime,5) .."/s. Hash:")
	else
		io.stdout:write("Download skipped. Hash:")
	end
	
	if checkSha1 then -- sha1 hash present?
		if item.file.hashes.sha1Hash then
			-- dont write hash twice
			_ = dltime>-1 and writeSha1Hash(path, dlpath, item.file.hashes.sha1Hash:lower())
			
			if genSha1(dlpath) == item.file.hashes.sha1Hash:lower() then
				print(" OK!")
			else
				print(" NOT OK!")
				return false
			end
		elseif item.file.hashes.quickXorHash then
			print(" NOCHK, quickXorHash: \"".. item.file.hashes.quickXorHash .."\"")
			_ = dltime>-1 and writeQuickXorHash(path, dlpath, item.file.hashes.quickXorHash)
		else
			print(" NOCHK, no remote hash or unknown")
		end
	end
	
	return true
end

function dlFolder(apiUrl, path, itemsDone)
	path = path or ("./" .. getFolderName(apiUrl) .. "/")
	path = sanitisePathLength(path, true)
	
	mkdir(path)
	
	itemsDone = itemsDone or 0
	
	local rootJson = wget(apiUrl)
	
	local tbl = JSON:decode(rootJson)
	assert(tbl, "No valid response received!")
	
	local itemsTotal = tbl["@odata.count"] or (#tbl.value ~= 0 and -#tbl.value or -1) --Sharepoint doesnt provide this
	print("Items total: ", itemsTotal)
	for n,item in pairs(tbl.value) do
		itemsDone = itemsDone + 1
		
		if sharepointDownload then checkRenewSharepointToken() end
		
		if item.folder then
	
			local ipath = path .. sanitisePathLength(item.name, true) .. "/"
			print("Folder: ".. item.name)
			mkdir(ipath)
			
			dlFolder(toShareFolderUrl(item.webUrl), ipath)

		elseif item.file then
			local dlpath = path .. sanitisePathLength(item.name)
			
			if skipDl == true then
				print("Skipping: ".. dlpath)
			else
				print("\n("..itemsDone.."/"..itemsTotal ..") File: ".. formatSize(item.size, 5) .." ".. item.name)
				local tries, wgetRetries = 0, sharepointDownload and 1 or nil
				repeat
					tries = tries + 1
					local ok = apiDownloadFile(item, path, dlpath, wgetRetries)
					
					if ok then
						break
					elseif not ok and sharepointDownload then
						sleepRetry(tries)
						checkRenewSharepointToken()
					end
				until tries == 15
				if tries == 15 then
					print("Giving up on file: ".. dlpath)
				end
			end
		else
			print("Unknown item type: ".. item.name)
		end
	end

	if tbl["@odata.nextLink"] then
		print(">> ".. itemsDone .."/".. itemsTotal .." items processed. Fetching more...")
		dlFolder(tbl["@odata.nextLink"], path, itemsDone)
	end
end

function writeSha1Hash(where, filePath, hash)
	local f = io.open(where .. "sha1.exf", "a+")
	f:write(hash," ?SHA1*", filePath, "\n")
	f:close()
end

function writeQuickXorHash(where, filePath, hash)
	local f = io.open(where .. "quickxor.exf", "a+")
	f:write(hash," ?QUICKXOR*", filePath, "\n")
	f:close()
end

function checkRenewSharepointToken()
	if os.time() - sharepointTokenSince > 3000 then --max 3600s
		sharepointTokenSince = os.time()
		print("Sharepoint token about to expire, renewing!")
		wgetEx(url, "-q &>/dev/null")
	end
end

function encOnedriveUrl(url)
	local ub64 = to_base64(url)
	
	return "u!" .. ub64:gsub("=+$", ""):gsub("/", "_"):gsub("%+", "-")
end

function toShareUrl(webUrl)
	local encUrl = encOnedriveUrl(webUrl)
	return "https://api.onedrive.com/v1.0/shares/".. encUrl .."/root/"
end

function toShareFolderUrl(webUrl)
	return toShareUrl(webUrl) .. "children"
end

url = arg[1]
if not url then
	error("Share URL not provided as first arg!")
else
	url = url:gsub("^['\"]", ""):gsub("['\"]$", "")
end

sharepointDownload = false
sharepointTokenSince = 0
extractLinksOnly = false
if arg[2] then
	local a = arg[2]
	if a == "-p" then
		print("Provided URL: ".. url)
		print("Share URL: ".. toShareUrl(url))
		os.exit(0)
	elseif a == "-s" then
		print("Forcing Sharepoint-downloader!")
	elseif a == "-w" then
		extractLinksOnly = true
		print("Printing wget download commands, no download!")
		print("The current launch options are:")
		print("wget ".. wgetOptions)
		print("\n")
	end
end

print("Downloading from this share url: ".. url)

if sharepointDownload or url:find("sharepoint.com", 1, true) then
	if not sharepointDownload then
		print("Auto-detected as Sharepoint-URL, good")
		sharepointDownload = true
	end
	
	-- original shortened URL sets you the access token
	-- the redirection URL doesnt and would ask for login (when tried in browser)
	local redirectUrl = getRedirect(url)
	if not redirectUrl then
		error("Couldn't fetch 302 Redirect URL!")
	end
	if redirectUrl:match("Throttle%.htm#?%d-$") then
		error("Hit a rate-limit, we are throttled!")
	end
	sharepointTokenSince = os.time()
	
	--print("\nSharepoint-302: ".. redirectUrl)
	
	local idFolder = redirectUrl:match("id=([^&]+)")
	print("Sharepoint base folder:", (idFolder:gsub("%%2[Ff]", "/")))
	
	-- https://<someSite>.sharepoint.com/:f:/g/personal/<DRIVENAME>/<HashStr>
	local driveName = url:match(".+(/.-/)") -- strip last URL part
	local driveUrl = redirectUrl:match(".-".. driveName)
	local rootFolderUrlArg, rootFolderPath = redirectUrl:match("([&?]RootFolder=([^&]+))")
	print("Sharepoint assuming drive URL:", driveUrl, "\n")
	if rootFolderUrlArg then
		print("Sharepoint RootFolder: ".. rootFolderPath:gsub("%%2[Ff]", "/"))
	else
		print("Sharepoint: Constructing RootFolder param from URL data (meh)")
		rootFolderPath = idFolder:gsub("%%2[Ff]$","") -- idFolder needs to be stripped of leading and trailing slashes
		rootFolderUrlArg = "&RootFolder="..rootFolderPath
	end
	
	local restListUrl = driveUrl .. "_api/web/GetListUsingPath(DecodedUrl=@a1)/RenderListDataAsStream?@a1=%27" .. idFolder .. "%27&TryNewExperienceSingle=TRUE" .. (rootFolderUrlArg or "")
	
	-- Get full metadata, RenderOptions is a bitmask
	-- https://github.com/SharePoint/PnP-JS-Core/pull/742
	-- https://docs.microsoft.com/en-us/dotnet/api/microsoft.sharepoint.client.renderlistdataoptions
	-- RenderOptions:
	-- 464647 - default: [listingJson.ListData.Row]
	-- 464386 - custom, less junk: [listingJson.Row]
	local listingPage = wgetEx(restListUrl, 
		"-O - --header='Content-Type: application/json;odata=verbose' "
		.. "--header='Accept: application/json;odata=verbose' "
		.. [[ --post-data='{"parameters":{"__metadata":{"type":"SP.RenderListDataParameters"},"RenderOptions":464386,"AllowMultipleValueFilterForTaxonomyFields":true,"AddRequiredFields":true}}']])
	local listingJson = JSON:decode(listingPage)
	assert(listingJson, "Couldnt fetch/parse Sharepoint listing!")
	
	-- iterate over individual items/folders
	for n, item in pairs(listingJson.Row) do
		-- eg "https://<someSite>.sharepoint.com/_api/v2.0/drives/b!<LongDriveId>/items/<ItemId>?version=Published"
		local itemApiUrl = (item[".spItemUrl"]:gsub("%?.+", "")):gsub("\\u002f", "/") -- remove URL params
		
		if item.FSObjType == "0" then
			-- file
			-- genius me. The previous version DID NOTHING if it encountered a file. -> download that file? nah... silently do NOTHING. And this entire script was supposed to be a downloader. I'm great.
			
			local folderName = item.FileRef:match(".+/(.+)/") -- don't sort by Editor, not intended here.
			-- https://<site>.sharepoint.com/_api/v2.0/drives/b!<drive-b64>/items/<item-upper-b64>?version=Published
			local parentFolderApiUrlChildren = listingJson.CurrentFolderSpItemUrl:gsub("([^?]+)(.*)$","%1/children%2") -- appends /children before URL params
			print("Downloading the entire folder you have pointed to (via URL)")
			dlFolder(parentFolderApiUrlChildren, folderName .."/")
			
			print("!!!")
			print("Please check that all files and folders were downloaded recursively. I lack testing to know for sure, but all should be fine")
			print("Report back here: https://gitlab.com/VADemon/onedeath/-/issues/7")
			print("Assuming everything was already downloaded, quitting...")
			print("!!!")
			break
		elseif item.FSObjType == "1" then
			-- folder
			local editorName = item.Editor and item.Editor[1] and item.Editor[1].title
								or ((item.FileRef:match("/?(.-/[^/]+)"):gsub("[/ ]", "")) or "unknown")
			local folderName = item.FileLeafRef
			if folderName:find("2001") then
				dlFolder(itemApiUrl .."/children", editorName .."/".. folderName .. "/")
			end
		elseif item.FSObjType == "-1" then
			print("Encountered FSObjType = -1! This is an invalid object! UniqueID=".. (item.UniqueId or "N/A") .. ". Skipping...")
		elseif item.FSObjType == "2" then
			print("Encountered FSObjType = 2! This is a \"web site\". UniqueID=".. (item.UniqueId or "N/A") .. ". Skipping...")
		end
		--print(n, item[".spItemUrl"])
	end
else
	dlFolder(toShareFolderUrl(url))
end
