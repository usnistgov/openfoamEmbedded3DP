import foldermover as fm
import folderparser as fp
import sys
LOGGERDEFINED = False
LOGGERDEFINED = fp.openLog('folderparser.log', LOGGERDEFINED)

fm.doneFolder(sys.argv[1], 2.5)