
from scr.read_data import read_config, read_data


def main():
    cfg = read_config('config.yml')
    data = read_data(cfg)

if __name__ == "__main__":
    main()
