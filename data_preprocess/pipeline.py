import download_glycosylation
import download_s_nitrosylation
import merge


def main():
    download_glycosylation.main()
    download_s_nitrosylation.main()
    merge.main()

    print("\nEverything worked! :)\n")


if __name__ == '__main__':
    main()
