import cv2

def count_frames(video_path):
    cap = cv2.VideoCapture(video_path)
    if not cap.isOpened():
        raise IOError(f"Cannot open {video_path}")
    # CAP_PROP_FRAME_COUNT returns total number of frames
    total = int(cap.get(cv2.CAP_PROP_FRAME_COUNT))
    cap.release()
    return total

if __name__ == "__main__":
    path = path = r"C:\Users\HP\OneDrive - NTNU\Desktop\Master\Code\PostProcess\Resultater\papirtest6 - 21 sek.avi"

    print(f"Total frames: {count_frames(path)}")

