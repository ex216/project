import random

def main():
    mag = [1,2,3,4,5,6]
    while True:
        input("Press \"Enter\" to pull the trigger")
        i = random.choice(mag)
        if i == 6:
            print("You died")
            break
        elif i in mag:
            mag.remove(i)
            print(mag)

if __name__ == "__main__":
    main()
