import sys
import getopt

def main(argv):
    ''' main function of the converter '''
    inputname = ''
    outputname = ''

    # Parsing input file name
    try:
        opts, args = getopt.getopt(argv,"hi:",["ifile="])
    except getopt.GetoptError:
        print('OFFtoNGB.py -i <inputfile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('OFFtoNGB.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            inputname = arg

    inputfile = open(inputname, 'r')
    outputname = inputname.replace("models", "data").replace("off", "data")
    outputfile = open(outputname, 'w')

    inputfile.readline()
    n_vert, _, _ = list(map(int, inputfile.readline().split()))
    output_content = str(n_vert) + "\n"
    for _ in range(n_vert):
        output_content = output_content + inputfile.readline()
    inputfile.close()
    outputfile.write(output_content)
    outputfile.close()

    print("The data file has been created in " + outputname + ".")

if __name__ == "__main__":
    main(sys.argv[1:])
