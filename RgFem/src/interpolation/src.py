import os

def rename(base):
    for root, ds, fs in os.walk(base):
        for f in fs:
            if f.endswith('.C'):
                fullname = os.path.join(root, f)
                newName = fullname[:-2] + '.cpp'
                os.rename(fullname,newName)

def findAllFile(base, end):
    for root, ds, fs in os.walk(base):
        for f in fs:
            if f.endswith(end):                
                yield f

def main():
    base = './'
    rename(base)
    out = './src.txt'
    fout = open(out,"w")
    for i in findAllFile(base, '.h'):
        fout.write(i + '\n')
    for i in findAllFile(base, '.cpp'):
        fout.write(i + '\n')
    fout.close()

if __name__ == '__main__':
    main()