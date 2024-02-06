from ilmn.pelops import notifications


class TestNotification:
    def test_notify(self):
        service = notifications.SimpleNotificationService()
        service.notify("Hello")
        assert service.get_notifications() == ["Hello"]
