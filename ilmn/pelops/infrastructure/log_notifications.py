import functools
import logging

from ilmn.pelops import notifications


class LogNotificationService(notifications.NotificationService):
    def notify(self, message: str) -> None:
        logging.info(message)


class LogNotificationServiceFactory(notifications.NotificationServiceFactory):
    @functools.lru_cache(maxsize=1)
    def build(self) -> notifications.NotificationService:
        logging.basicConfig(level=logging.INFO)
        result = LogNotificationService()
        return result
