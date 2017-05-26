##################
### CHANNEL TYPE CHECKS
###################	
def elecChannelsList():
	elecChannels = ["Ex", "Ey"]
	return elecChannels

def isElectric(chan):
	if chan in elecChannelsList():
		return True
	return False

def magChannelsList():
	magChannels = ["Hx", "Hy", "Hz", "Bx", "By", "Bz"]
	return magChannels

def isMagnetic(chan):
	if chan in magChannelsList():
		return True
	return False

def consistentChans(chan):
	standardChans = ["Hx", "Hy", "Hz", "Ex", "Ey"]
	if chan in standardChans:
		return chan
	if chan == "Bx":
		return "Hx"
	if chan == "By":
		return "Hy"
	if chan == "Bz":
		return "Hz"
	# otherwise return chan
	return chan