from biaxial_analysis import run
import pickle
import sys

if __name__ == "__main__":
    if len(sys.argv) > 1:
        # Pick serialized object
        try:
            with open(sys.argv[1], 'rb') as f:
                obj = pickle.load(f)
        except IOError:
            raise Exception("File not found")
    else:
        raise Exception("No file provided")

    # Run
    run(obj)