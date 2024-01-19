import os
from datetime import datetime

class StorageUtils:
    def __init__(self,
                storage_dir:str="SEWAGE_Workspace",
                time_stamp:bool=False,
                ):
        self.storage_dir = os.path.realpath(storage_dir)
        self.time_stamp = time_stamp
        self.date_time = datetime.now().strftime("%Y%m%d_%H%M%S")

    def get_storage_pathway(self):
        return self.storage_dir
    
    def add_time_stamp(self):
            if self.time_stamp:
                self.storage_dir = f"{self.storage_dir}_{self.date_time}"
            else:
                self.storage_dir = f"{self.storage_dir}"

    def create_storage_dir(self):
        try:
            os.mkdir(self.storage_dir)
        except FileExistsError:
            raise FileExistsError(f"Directory '{self.storage_dir}' already exists.")
        except FileNotFoundError:
            raise FileNotFoundError()
        except Exception as e:
            print(f"An error occurred: {e}")
            raise e