timePerPair = (13 / 6429)
totalPairs = (505 - 31 + 664 - 518) ** 2
ts = totalPairs * timePerPair
print("total secs: " + str(ts))
print("total mins: " + str(ts / 60))
print("total hours: " + str(ts / 3600))
print("total days: " + str(ts / (3600 * 24)))
