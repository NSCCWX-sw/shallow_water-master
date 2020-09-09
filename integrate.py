import os
import sys

def main():
    data = []
    temp = []
    files = []
    rows = int(sys.argv[1])
    cols = int(sys.argv[2])
    count = 0
    row = 0
    # aggregate files
    for i in range(rows*cols):
        file_name = 'eta01000_'+ str(i) + '.dat'
        files.append(file_name)
    # check validity
    for each_file in files:
        path = os.path.join('eta', each_file)
        if os.path.exists(path) is False:
            print each_file, "doesn't exist.\n"

    for idx, each_file in enumerate(files, 0):
        path = os.path.join('eta', each_file)
        print(path)
        with open(path, 'r') as f1:
            list1 = f1.readlines()
            row = int(count/cols)
            if idx%cols==0:
                for i in range(0, len(list1)):
                    data.append(list1[i].rstrip('\n'))
            else:
                for i in range(0, len(list1)):
                    temp = data[i+len(list1)*row]
                    temp += ' '
                    temp += list1[i].rstrip('\n')
                    data[i+len(list1)*row] = temp
            count += 1
    
    print("It comes here")
    path = os.path.join('eta', "final_data.dat")
    f_write = open(path, "w+")
    for i in range(len(data)):
        f_write.write(data[i]+'\n')
    f_write.close()  

if __name__ == "__main__":
    main()  
