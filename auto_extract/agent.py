import anthropic
import time
from openai import AzureOpenAI, OpenAI
from pypdf import PdfReader
from tqdm import tqdm
from abc import ABC, abstractmethod
from typing import List
import numpy as np

from .config import Config

class BaseClient(ABC):
    def __init__(self, pdf_path: str = None):
        self.api_key = Config().API_KEY
        self.api_base_url = Config().API_BASE_URL
        self.model = Config().MODEL
        self.tool_model = Config().TOOL_MODEL
        self.messages = []
        
        self.pdf_path = pdf_path
        self.extracted_text = ''
        if self.pdf_path:
            if self.pdf_path.endswith('.txt'):
                with open(self.pdf_path, 'r') as f:
                    self.extracted_text = f.read()
            else:
                self.extracted_text = self.extract_text_from_pdf(self.pdf_path)
    
    def extract_text_from_pdf(self, pdf_path: str) -> str:
        reader = PdfReader(pdf_path)
        text = ''
        for page in tqdm(reader.pages, desc='Extracting text from PDF'):
            text += page.extract_text()
            text += '\n'
        return text
    
    def update_messages(func):
        def wrapper(self, *args, **kwargs):
            resp = func(self, *args, **kwargs)
            self.messages.append({"role": "assistant", "content": resp})
            return resp
        return wrapper
    
    @update_messages
    def initiate_propmt(self, 
                        initial_message: str = Config().get_prompt('SYSTEM_PROMPT'),
                        user_message: str = Config().get_prompt('USER_ARTICLE_PROMPT'),
                        ):
        
        if self.extracted_text == '':
            raise ValueError('No PDF file provided')

        self.messages.append({"role": "system", "content": initial_message})
        user_message += self.extracted_text
        self.messages.append({"role": "user", "content": user_message})
        
        resp = self._retrieve(self.messages)
        
        return resp
    
    @update_messages
    def chat(self, new_message: str, max_tokens = 1000):

        self.messages.append({"role": "user", "content": new_message})
        resp = self._retrieve(self.messages, max_tokens=max_tokens)
        
        return resp
    
    def clear_intermediate_messages(self):
        # only keep the system message, extract text message, response message
        if len(self.messages) > 3:
            self.messages = self.messages[:3]
        else:
            pass
    
    @abstractmethod
    def create_client(self):
        pass
    
    @abstractmethod
    def _retrieve(self):
        pass
    
    @abstractmethod
    def _tool_retrieve(self):
        pass

class Openai(BaseClient):
    def __init__(self, pdf_path: str = None):
        super().__init__(pdf_path)
        self.create_client()

    def create_client(self):
        self.client = OpenAI(api_key=self.api_key, base_url=self.api_base_url)

    # add retrieve method   
    def _retrieve(self,
                 messages: list, 
                 max_tokens: int = 1000,
                 max_retries: int = 3
                 ):

        retry_times = 0
        while retry_times < max_retries:
            try:
                completion = self.client.chat.completions.create(
                    model=self.model,
                    messages=messages,
                    max_tokens=max_tokens
                )
                break
                
            except Exception as e:
                retry_times += 1
                time.sleep(5*retry_times)
                print(f"Error: {e}, retrying {retry_times} times")
                if retry_times == max_retries:
                    raise e
                continue
        
        return completion.choices[0].message.content
    
    def _tool_retrieve(self,
                       messages: list,
                       max_tokens: int = 1000,
                       max_retries: int = 3
                       ):
        
        retry_times = 0
        while retry_times < max_retries:
            try:
                completion = self.client.chat.completions.create(
                    model=self.tool_model,
                    messages=messages,
                    max_tokens=max_tokens
                )
                break
                
            except Exception as e:
                time.sleep(5*retry_times)
                retry_times += 1
                print(f"Error: {e}, retrying {retry_times} times")
                if retry_times == max_retries:
                    raise e
                continue
            
        return completion.choices[0].message.content
    
class Claude3(BaseClient):
    def __init__(self, pdf_path: str = None):
        super().__init__(pdf_path)
        self.create_client()
        
    def create_client(self):
        self.client = anthropic.Anthropic(api_key=self.api_key)
    
    # add retrieve method
    def _retrieve(self, 
                 messages: list, 
                 max_tokens: int = 1000,
                 max_retries: int = 3
                 ):
        
        retry_times = 0
        while retry_times < max_retries:
            try:
                completion = self.client.messages.create(
                    model=self.model,
                    messages=messages,
                    max_tokens=max_tokens
                )
                break
                
            except Exception as e:
                time.sleep(5*retry_times)
                retry_times += 1
                print(f"Error: {e}, retrying {retry_times} times")
                if retry_times == max_retries:
                    raise e
                continue
        
        return completion.content[0].text
    
    def _tool_retrieve(self,
                        messages: list,
                        max_tokens: int = 1000,
                        max_retries: int = 3
                        ):
          
        retry_times = 0
        while retry_times < max_retries:
            try:
                completion = self.client.messages.create(
                    model=self.tool_model,
                    messages=messages,
                    max_tokens=max_tokens
                )
                break
                
            except Exception as e:
                retry_times += 1
                time.sleep(5*retry_times)
                print(f"Error: {e}, retrying {retry_times} times")
                if retry_times == max_retries:
                    raise e
                continue
        
        return completion.content[0].text

    def update_messages(func):
        def wrapper(self, *args, **kwargs):
            resp = func(self, *args, **kwargs)
            self.messages.append({"role": "assistant", "content": resp})
            return resp
        return wrapper

    @update_messages
    def initiate_propmt(self, 
                        initial_message: str = Config().get_prompt('SYSTEM_PROMPT'),
                        user_message: str = Config().get_prompt('USER_ARTICLE_PROMPT'),
                        ):
        
        if self.extracted_text == '':
            raise ValueError('No PDF file provided')

         # claude3 does not have a system message
        initial_message += user_message
        initial_message += self.extracted_text
        self.messages.append({"role": "user", "content": initial_message})
        
        resp = self._retrieve(self.messages)
        
        return resp
    
def get_cell_type_embedding_by_llm(cell_types: List[str]) -> List[np.ndarray]:
    """
    Get cell type embeddings by using the OpenAI API.
    """
    embedding_api_key = Config().EMBEDDING_API_KEY
    azure_endpoint = Config().EMBEDDING_ENDPOINT
    client = AzureOpenAI(
    api_key=embedding_api_key, api_version="2024-02-01", azure_endpoint=azure_endpoint
    )
    response = client.embeddings.create(
        model=Config().EMBEDDING_MODEL, input=cell_types
    )
    emb = [x.embedding for x in response.data]
    return emb