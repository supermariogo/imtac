import pandas as pd
from openpyxl import load_workbook
import os
import math
import statistics
from numpy import mean, std
import sys


hash = dict() # key is 
final = pd.DataFrame()


def firstOccurance(a, x):
    lo = 0
    hi = len(a) - 1

    while lo < hi - 1:
        mid = (lo + hi) // 2

        if a[mid] < x:
            lo = mid
        elif a[mid] > x:
            hi = mid
        else:
            hi = mid
    if a[hi] == x:
    	return hi + 1
    elif a[lo] == x:
    	return lo + 1
    return -1


def getAbundanceColumnName(df, sheetName):
	columnNames = list(df.columns.values)
	for columnName in columnNames:
		if 'Abundance Ratio' in columnName and sheetName.split('-')[0] in columnName and sheetName.split('-')[1] in columnName:
			return columnName

def processFile(fileName, sheetName, carriedFiled):
	df = pd.read_excel(fileName, dtype=str)
	df = df[df["Protein FDR Confidence: Combined"] == 'High']
	orignalColumns = ['Accession', 'Gene Symbol', '# Unique Peptides']


	# result_df = pd.DataFrame(columns = ['Accession', 'Gene Symbol', '# Unique Peptides', 'Abundance', 'Tag', 'IMTAC ID', 'Date', 'Concentration', 'Cell Line', 'Probe','Raio N', 'Score', 'Rank', 'Log2', 'Median'])
	result_df = pd.DataFrame(columns = orignalColumns)
	for orignalColumn in orignalColumns:
		result_df[orignalColumn] = df[orignalColumn]
	# abundance Ratio
	abundanceRatioName = getAbundanceColumnName(df, sheetName)
	result_df["Abundance Ratio"] = df[abundanceRatioName]
	abundanceRatioName = "Abundance Ratio"
	# carried values
	for key, value in carriedFiled.items():
		result_df[key] = value
	#log2
	result_df['Log2'] = result_df.apply(lambda row: math.log2(float(row[abundanceRatioName])), axis = 1)

	#Set log2 to nan if gene Symbol start with KRT
	for i, j in result_df.iterrows(): 
		if j['Gene Symbol'].startswith('KRT'):
			j['Log2'] = 'nan'

	# score
	log2Floats = [x for x in list(result_df['Log2']) if str(x) != 'nan']
	meanValue = mean(log2Floats)
	stdValue = statistics.stdev(log2Floats)
	#print("std: " + str(stdValue) + "mean: " + str(meanValue))
	result_df['Score'] = result_df.apply(lambda row: (row['Log2'] - meanValue)/stdValue, axis = 1)



	#rank
	allScores = sorted([x for x in list(result_df['Score']) if str(x) != 'nan'])
	result_df['Rank'] = result_df.apply(lambda row: firstOccurance(allScores, row['Score']), axis = 1)


	#Calculate Median
	all_ratios = []
	for ratio in result_df[abundanceRatioName]:
		try:
			if str(ratio) == 'nan':
				continue
			all_ratios.append(float(ratio))
		except:
			print("AbundanceRatio is not convertible")
	median = statistics.median(all_ratios)
	result_df['Median'] = median

	# Raion N
	result_df['Ratio N'] = result_df.apply(lambda row: float(row[abundanceRatioName])/median, axis = 1)

	#Position
	result_df['Position'] = result_df.apply(lambda row: str(row['Rank']) + "/" + str(len(all_ratios)), axis = 1)

	#append to final
	global final
	final = final.append(result_df)[result_df.columns.tolist()]
	return result_df

def writeToFile(fileName, sheetName, result_df):
	# https://stackoverflow.com/questions/40385689/add-a-new-sheet-to-a-existing-workbook-in-python
	# book = load_workbook("/Users/zemao/Desktop/result/"+fileName)
	# writer = pd.ExcelWriter("/Users/zemao/Desktop/result/"+fileName, engine='openpyxl') 
	# writer.book = book
	# writer.sheets = dict((ws.title, ws) for ws in book.worksheets)
	if not os.path.exists('processed'):
		os.makedirs("processed")
	df = pd.read_excel(fileName, header=None)
	writer = pd.ExcelWriter("processed/" + fileName, engine = 'openpyxl')
	df.to_excel(writer, "Proteins", header=None, index=None)
	result_df.to_excel(writer, sheetName, index=None)
	writer.save()
	writer.close()

def writeAllTofile(path):
	writer = pd.ExcelWriter(path, engine = 'openpyxl')
	final.to_excel(writer, "all")
	writer.save()
	writer.close()



if __name__ == '__main__':
	try:
		inputFile = "log.csv" if len(sys.argv) <= 1 else sys.argv[1]
		df = pd.read_csv(inputFile, header=None, dtype=str)
		fileNames = os.listdir('.')
		for i in range(len(df)) : 
			#print(list(df.loc[i]))
			s = df.loc[i, 1].replace("S", "")
			if len(s.split('-')) != 2:
				print("Skip", df.loc[i, 1])
				continue

			fileNameKeyword = s.split('-')[0] + "-S" + s.split('-')[1]
			number = fileNameKeyword.split('-')[0]
			sn = fileNameKeyword.split('-')[1]
			people = df.loc[i, 2]
			
			small = str(int(df.loc[i, 3])-1)
			big = df.loc[i, 3]
			sheetName = small+"-"+big
			print("Process " + number + " " +sn + " " + people)
			carriedFiled = dict()
			carriedFiled['Tag'] = big +'/' + small
			carriedFiled["IMTAC ID"] = "IMTAC%s-%s-%s"%(number, sn, people)
			carriedFiled["Date"] = df.loc[i, 0].split(" ")[0]
			carriedFiled["Concentration"] = str(df.loc[i, 6])+"/"+str(df.loc[i, 8])
			carriedFiled["Cell Line"] = df.loc[i, 4]
			carriedFiled["Probe"] = df.loc[i, 7]

			for fileName in fileNames:
				if number in fileName and sn in fileName and people in fileName and "~" not in fileName:
					print("Working on file " + fileName)
					result_df = processFile(fileName, sheetName, carriedFiled)
					print("Writing to file " + fileName)
					writeToFile(fileName, sheetName, result_df)
			print("------------")
		print("All done, generating all.")
		writeAllTofile("processed/all.xlsx")
	except Exception as e:
		print(e)
	input('Program done. Enter anything to close this window.')
